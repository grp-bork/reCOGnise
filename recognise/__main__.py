import argparse
import json
import os
import pathlib
import re
import subprocess
import sys

from collections import Counter

from .db.queries import get_sequences_from_cluster


SPECI_COGS = """
	COG0012 COG0016 COG0018 COG0048 COG0049 COG0052 COG0080 COG0081
	COG0085 COG0087 COG0088 COG0090 COG0091 COG0092 COG0093 COG0094
	COG0096 COG0097 COG0098 COG0099 COG0100 COG0102 COG0103 COG0124
	COG0172 COG0184 COG0185 COG0186 COG0197 COG0200 COG0201 COG0202
	COG0215 COG0256 COG0495 COG0522 COG0525 COG0533 COG0541 COG0552
""".strip().split(" ")

MOTU_COGS = set(
    """
	COG0012 COG0016 COG0018 COG0172 COG0215
	COG0495 COG0525 COG0533 COG0541 COG0552
	""".strip().split(" ")
)

COGS = {
	cog: cog in MOTU_COGS
	for cog in SPECI_COGS
}


def main():
	
	ap = argparse.ArgumentParser()
	ap.add_argument("genome_id", type=str)
	ap.add_argument("--genes", type=str)
	ap.add_argument("--proteins", type=str)
	ap.add_argument("--genome", type=str)
	ap.add_argument("cog_db", type=str)
	ap.add_argument("--cpus", type=int, default=4)
	ap.add_argument("--output_dir", "-o", type=str, default="recognise_out")
	ap.add_argument("--dbcred", type=str)
	
	args = ap.parse_args()

	genome_present, genes_present, proteins_present = (f is not None and os.path.isfile(f) for f in (args.genome, args.genes, args.proteins))
	genes, proteins = None, None
	
	pathlib.Path(args.output_dir).mkdir(exist_ok=True, parents=True)
	
	if genome_present:
		if genes_present or proteins_present:
			raise ValueError("Please specify either a genome or a gene/protein set combination.")
		
		# call prodigal
		proteins = os.path.join(args.output_dir, f"{args.genome_id}.faa")
		genes = os.path.join(args.output_dir, f"{args.genome_id}.ffn")
		prodigal_proc = subprocess.run(
			[
				# prodigal -i \$(basename ${genome_fna} .gz) -f gff -o ${genome_id}/${genome_id}.gff -a ${genome_id}/${genome_id}.faa -d ${genome_id}/${genome_id}.ffn
				"prodigal",
				"-i",
				args.genome,
				"-a",
				proteins,
				"-d",
				genes,
			], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
		)
		# out, err = prodigal_proc.communicate()

		if prodigal_proc.returncode != 0:
			raise ValueError(f"<pre>Prodigal error\n\n{prodigal_proc.stdout}</pre>")


	elif genes_present and proteins_present:
		genes, proteins = args.genes, args.proteins
	elif genes_present:
		raise ValueError("Missing protein set, please specify with --proteins.")
	elif proteins_present:
		raise ValueError("Missing gene set, please specify with --genes.")

	try:
		dbstr = json.load(open(args.dbcred, "rt")).get("DB_STR")
	except:
		dbstr = None

	# cog_dir = os.path.join(args.output_dir, "cogs")
	# cog_dir = f"{args.genome_id}_cogs"
	pathlib.Path(cog_dir).mkdir(exist_ok=True, parents=True)

	# fetchMG.pl -o \${bin_id}_cogs -t 5 -m extraction -d genecalls/\${bin_id}.extracted.fna genecalls/\${bin_id}.extracted.faa
	fetchmg_cmd = [
		"fetchMGs.pl",
		"-o", cog_dir,
		"-t", f"{args.cpus}",
		"-m", "extraction",
		"-x", "/usr/bin",
		"-d", f"{genes}",
		f"{proteins}",
	]
	fetchmg_proc = subprocess.run(
		fetchmg_cmd,
		stdout=subprocess.PIPE, stderr=subprocess.PIPE,     
	)    

	# out, err = fetchmg_proc.communicate()
	if fetchmg_proc.returncode != 0:
		raise ValueError(f"<pre>fetchMGs error\n\n{' '.join(fetchmg_cmd)}\n\n{fetchmg_proc.stdout.decode()}</pre>")


	# print(out.decode())

	align_file = os.path.join(cog_dir, "dummy.fna")

	speci_cog_d = {}
	speci_header = None

	for cog in COGS:
		cog_file = os.path.join(cog_dir, f"{cog}.fna")
		if os.path.isfile(cog_file):
			with open(align_file, "wt") as aln_file, open(cog_file, "rt") as cog_in:
				for line in cog_in:
					if line[0] == ">":
						line = re.sub(r"$", f"  # {cog} {args.genome_id}", line.strip())
					print(line.strip(), file=aln_file)


			mapseq_pr = subprocess.Popen(
				[
					"mapseq",
					align_file,
					os.path.join(args.cog_db, f"{cog}.fna"),
					os.path.join(args.cog_db, f"{cog}.specI.tax"),					
				],
				stdout=subprocess.PIPE, stderr=subprocess.PIPE,
			)

			out, err = mapseq_pr.communicate()

			out = out.decode().strip().split("\n")
			if speci_header is None:				
				speci_header = [line.strip().split("\t") for line in out if line[0] == "#"]
				speci_header[-1].insert(0, "cog")

			speci_cog_d[cog] = [line.strip().split("\t") for line in out if line[0] != "#"]
			for line in speci_cog_d[cog]:
				line.insert(0, cog)
			
			# cat mapseq/\${bin_id}/speci/* | sed '3,\${ /^#/d }' > ${sample_id}/\${bin_id}.speci.assignments

			# break

	cogs_out = open(os.path.join(args.output_dir, f"{args.genome_id}.cogs.txt"), "wt")
	speci_out = open(os.path.join(args.output_dir, f"{args.genome_id}.specI.txt"), "wt")
	speci_status_out = open(os.path.join(args.output_dir, f"{args.genome_id}.specI.status"), "wt")
	seqfile = os.path.join(args.output_dir, f"{args.genome_id}.specI.ffn.gz")
	open(seqfile, "wt").close()

	with cogs_out, speci_out, speci_status_out:
		specis = Counter()
		if speci_header:
			print(*("\t".join(line) for line in speci_header), sep="\n", file=cogs_out)
		for cog in COGS:
			cog_lines = speci_cog_d.get(cog)
			if cog_lines is not None:
				for line in cog_lines:
					print("\t".join(line), file=cogs_out)
					specis[line[14]] += 1
		
		
		speci_counts = specis.most_common()
		print(speci_counts)
		if not speci_counts:
			print("Warning: could not find any markers. Aborting.")
			print("NO_MARKERS", file=speci_status_out)
		else:
			first, second = (speci_counts[:2]) if len(speci_counts) > 2 else (speci_counts[0], None)

			if (second is None) or (first[1] > second[1]):
				n_seqs = -1
				if dbstr is not None:					
					n_seqs = get_sequences_from_cluster(dbstr, first[0], seqfile)

				print(first[0], file=speci_out)

				if not n_seqs:
					print("Warning: specI cluster is too small. Aborting.")
					print("SPECI_SIZE_INSUFFICIENT", file=speci_status_out)
				else:
					print("OK", file=speci_status_out)
					open(speci_status_out.name + ".OK", "wt").close()

			else:
				print(f"Warning: cannot determine consensus specI. first={first} second={second}. Aborting.")
				print("NO_CONSENSUS", file=speci_status_out)






	"""
	for cog in \${specicogs[@]}; do
			if [[ -s \${bin_id}_cogs/\${cog}.fna ]]
			then
				sed -i '/>/ s/\$/ '" # \${cog} \${bin_id}"'/' \${bin_id}_cogs/\${cog}.fna
				((i=++i%6)) || wait
				mapseq \${bin_id}_cogs/\${cog}.fna \${cogdir}\${cog}.fna \${cogdir}\${cog}.specI.tax | sed "s/#query/#cog\tquery/" | sed "1,2 ! s/^/\${cog}\t&/" > mapseq/\${bin_id}/speci/\${cog} &
			else
				touch mapseq/\${bin_id}/speci/\${cog}
			fi
		done
	"""
		
	
	

if __name__ == "__main__":
	main()

	"""
	mkdir ${sample_id}
	mkdir mapseq

	rsync -av ${params.cog_ref_dir} \$TMPDIR
	cogdir=\$TMPDIR/

	specicogs=(COG0012 COG0016 COG0018 COG0048 COG0049 COG0052 COG0080 COG0081
			   COG0085 COG0087 COG0088 COG0090 COG0091 COG0092 COG0093 COG0094
			   COG0096 COG0097 COG0098 COG0099 COG0100 COG0102 COG0103 COG0124
			   COG0172 COG0184 COG0185 COG0186 COG0197 COG0200 COG0201 COG0202
			   COG0215 COG0256 COG0495 COG0522 COG0525 COG0533 COG0541 COG0552)
	motu_cogs=(COG0012 COG0016 COG0018 COG0172 COG0215
			   COG0495 COG0525 COG0533 COG0541 COG0552)
	i=0

	echo -n "bin" > ${sample_id}/${sample_id}.cog_count
	for cog in \${specicogs[@]}; do
		echo -n "\t"\${cog} >> ${sample_id}/${sample_id}.cog_count
	done

	for bin in bins/*
	do
	  bin_id=\${bin:5:\${#bin}-8}
	  if [[ -s genecalls/\${bin_id}.extracted.fna  ]]
	  then
		fetchMG.pl -o \${bin_id}_cogs -t 5 -m extraction -d genecalls/\${bin_id}.extracted.fna genecalls/\${bin_id}.extracted.faa

		mkdir mapseq/\${bin_id}
		mkdir mapseq/\${bin_id}/motu
		mkdir mapseq/\${bin_id}/speci

		# SPECI
		for cog in \${specicogs[@]}; do
			if [[ -s \${bin_id}_cogs/\${cog}.fna ]]
			then
				sed -i '/>/ s/\$/ '" # \${cog} \${bin_id}"'/' \${bin_id}_cogs/\${cog}.fna
				((i=++i%6)) || wait
				mapseq \${bin_id}_cogs/\${cog}.fna \${cogdir}\${cog}.fna \${cogdir}\${cog}.specI.tax | sed "s/#query/#cog\tquery/" | sed "1,2 ! s/^/\${cog}\t&/" > mapseq/\${bin_id}/speci/\${cog} &
			else
				touch mapseq/\${bin_id}/speci/\${cog}
			fi
		done
		wait
		cat mapseq/\${bin_id}/speci/* | sed '3,\${ /^#/d }' > ${sample_id}/\${bin_id}.speci.assignments
		i=0
		# MOTU
		for cog in \${motu_cogs[@]}; do
			if [[ -s \${bin_id}_cogs/\${cog}.fna ]]
			then
				((i=++i%6)) || wait
				mapseq \${bin_id}_cogs/\${cog}.fna \${cogdir}\${cog}.mOTU.fna \${cogdir}\${cog}.mOTU.tax | sed "s/#query/#cog\tquery/" | sed "1,2 ! s/^/\${cog}\t&/" > mapseq/\${bin_id}/motu/\${cog} &
			else
				touch mapseq/\${bin_id}/motu/\${cog}
			fi
		done
		wait
		cat mapseq/\${bin_id}/motu/* | sed '3,\${ /^#/d }' > ${sample_id}/\${bin_id}.motu.assignments

		# SUMMARISE AND CLEANUP
		cat \${bin_id}_cogs/*.fna >> ${sample_id}/${sample_id}.marker_genes.fna
		taxonomy_annotation.py -m ${sample_id}/\${bin_id}.motu.assignments -s ${sample_id}/\${bin_id}.speci.assignments -o \${bin_id}.taxonomy_annotation -p \${bin_id}

		# GET COG COUNTS
		echo -n "\n"\${bin_id} >> ${sample_id}/${sample_id}.cog_count
		for cog in \${specicogs[@]}; do
			echo -n "\t"`grep -c \${cog} <(cat \${bin_id}_cogs/*.fna)` >> ${sample_id}/${sample_id}.cog_count
		done
		echo >> ${sample_id}/${sample_id}.cog_count
	  else
	  echo "No genecalls for \${bin_id}"
	  fi
	done

	echo "name\ttaxon_level\ttaxonomic_assignment\tassignment_confidence\tconfident_count\tother_confident_count\tunconfident_count\tother_unconfident_count\tunique_marker_genes\ttotal_marker_genes\tchimerism_score_all\tchimerism_score_confident\tari_all\tari_confident" > ${sample_id}.taxonomy_annotation.tsv
	cat *.taxonomy_annotation >> ${sample_id}.taxonomy_annotation.tsv

	gzip ${sample_id}/${sample_id}.marker_genes.fna
	mkdir ${sample_id}_assignments
	mv ${sample_id}/*.assignments ${sample_id}_assignments/
	tar -cvzf ${sample_id}/${sample_id}.assignments.tar.gz ${sample_id}_assignments/*.assignments
	"""
