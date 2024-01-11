import argparse
import itertools as it
import json
import logging
import multiprocessing as mp
import os
import pathlib
import re
import subprocess
import sys

from collections import Counter

from .db.queries import get_sequences_from_cluster
from .tasks.mapseq import task


logging.basicConfig(
	level=logging.DEBUG,
	format='[%(asctime)s] %(message)s'
)

logger = logging.getLogger(__name__)

SPECI_COGS = (
	"COG0012", "COG0016", "COG0018", "COG0048", "COG0049", "COG0052", "COG0080", "COG0081",
	"COG0085", "COG0087", "COG0088", "COG0090", "COG0091", "COG0092", "COG0093", "COG0094",
	"COG0096", "COG0097", "COG0098", "COG0099", "COG0100", "COG0102", "COG0103", "COG0124",
	"COG0172", "COG0184", "COG0185", "COG0186", "COG0197", "COG0200", "COG0201", "COG0202",
	"COG0215", "COG0256", "COG0495", "COG0522", "COG0525", "COG0533", "COG0541", "COG0552",
)

MOTU_COGS = {
	"COG0012", "COG0016", "COG0018", "COG0172", "COG0215",
	"COG0495", "COG0525", "COG0533", "COG0541", "COG0552",
}

COGS = {
	cog: cog in MOTU_COGS
	for cog in SPECI_COGS
}


def call_prodigal(genome, protein_file, gene_file):
	# prodigal -i \$(basename ${genome_fna} .gz) -f gff -o ${genome_id}/${genome_id}.gff -a ${genome_id}/${genome_id}.faa -d ${genome_id}/${genome_id}.ffn
	prodigal_proc = subprocess.run(
		[
			"prodigal",
			"-i",
			genome,
			"-a",
			protein_file,
			"-d",
			gene_file,
		],
		stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
	)
	
	if prodigal_proc.returncode != 0:
		raise ValueError(f"<pre>Prodigal error\n\n{prodigal_proc.stdout}</pre>")
	

def call_fetch_mgs(protein_file, gene_file, cog_dir, cpus):
	# fetchMG.pl -o \${bin_id}_cogs -t 5 -m extraction -d genecalls/\${bin_id}.extracted.fna genecalls/\${bin_id}.extracted.faa
	fetchmgs_cmd = [
		"fetchMGs.pl",
		"-o", cog_dir,
		"-t", f"{cpus}",
		"-m", "extraction",
		"-x", "/usr/bin",
		"-d", f"{gene_file}",
		f"{protein_file}",
	]		
	fetchmgs_proc = subprocess.run(
		fetchmgs_cmd,
		stdout=subprocess.PIPE, stderr=subprocess.STDOUT,     
	)
	    
	if fetchmgs_proc.returncode != 0:
		raise ValueError(f"<pre>fetchMGs error\n\n{' '.join(fetchmgs_cmd)}\n\n{fetchmgs_proc.stdout.decode()}</pre>")


# def call_mapseq(align_file, cog_db, cog, speci_header=None):
# 	mapseq_pr = subprocess.Popen(
# 		[
# 			"mapseq",
# 			align_file,
# 			os.path.join(cog_db, f"{cog}.fna"),
# 			os.path.join(cog_db, f"{cog}.specI.tax"),					
# 		],
# 		stdout=subprocess.PIPE, stderr=subprocess.PIPE,
# 	)

# 	out, err = mapseq_pr.communicate()

# 	msg, speci_header, speci_cog = None, None, None
# 	if mapseq_pr.returncode != 0:
# 		msg = err.decode().strip().split("\n")
# 	else:
# 		out = out.decode().strip().split("\n")

# 		if speci_header is None:				
# 			speci_header = [line.strip().split("\t") for line in out if line[0] == "#"]
# 			speci_header[-1].insert(0, "cog")

# 		speci_cog = [[cog] + (line.strip().split("\t")) for line in out if line[0] != "#"]	

# 	return msg, speci_header, speci_cog

# def task(cog_file, cog, genome_id, cog_db):
# 	with open(cog_file + ".align", "wt") as aln_file, open(cog_file, "rt") as cog_in:
# 		for line in cog_in:
# 			if line[0] == ">":
# 				line = re.sub(r"$", f"  # {cog} {genome_id}", line.strip())
# 			print(line.strip(), file=aln_file)
	
# 	msg, _, cog_lines = call_mapseq(cog_file + ".align", cog_db, cog)
# 	return msg, cog_lines


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
	ap.add_argument("--marker_set", type=str, choices=("full", "motus", "test"), default="motus")
	
	args = ap.parse_args()

	genome_present, genes_present, proteins_present = (f is not None and os.path.isfile(f) for f in (args.genome, args.genes, args.proteins))
	genes, proteins = None, None
	
	pathlib.Path(args.output_dir).mkdir(exist_ok=True, parents=True)
	
	if genome_present:
		logger.info("Running prodigal...")
		if genes_present or proteins_present:
			raise ValueError("Please specify either a genome or a gene/protein set combination.")
		
		proteins = os.path.join(args.output_dir, f"{args.genome_id}.faa")
		genes = os.path.join(args.output_dir, f"{args.genome_id}.ffn")

		call_prodigal(args.genome, proteins, genes)
		logger.info("prodigal finished.")

	elif genes_present and proteins_present:
		genes, proteins = args.genes, args.proteins
	elif genes_present:
		raise ValueError("Missing protein set, please specify with --proteins.")
	elif proteins_present:
		raise ValueError("Missing gene set, please specify with --genes.")

	try:
		# dbstr = json.load(open(args.dbcred, "rt")).get("DB_STR")
		db_d = json.load(open(args.dbcred, "rt")).get("progenomes3_db")
	except:
		db_d = {}

	user = db_d.get("username")
	host = db_d.get("host")
	pw = db_d.get("password")
	port = db_d.get("port")

	dbstr = f"mongodb://{user}:{pw}@{host}:{port}" if (user and host and pw and port) else None

	cog_dir = os.path.join(args.output_dir, "cogs")
	
	logger.info("Running fetchMGs...")
	call_fetch_mgs(proteins, genes, cog_dir, args.cpus)
	logger.info("fetchMGs finished.")

	align_file = os.path.join(cog_dir, "dummy.fna")

	speci_cog_d = {}
	speci_header = None
	specis = Counter()

	tasks = []
	for cog, is_motus_cog in COGS.items():
		if args.marker_set == "motus" and not is_motus_cog:
			continue
		if args.marker_set == "test" and len(tasks) == 3:
			break
		cog_file = os.path.join(cog_dir, f"{cog}.fna")
		if os.path.isfile(cog_file) and os.stat(cog_file).st_size:
			tasks.append((cog_file, cog, args.genome_id, args.cog_db, min(args.cpus, 4)))

	logger.info(f"Running {args.cpus // min(args.cpus, 4)} MAPseq processes on {len(tasks)} marker genes. marker_set={args.marker_set}...")

	with mp.Pool(args.cpus // min(args.cpus, 4)) as pool:
		# results = list(it.chain(*pool.starmap_async(task, tasks).get()))
		results = list(pool.starmap_async(task, tasks).get())

	logger.info("MAPseq finished.")

	print(results)

	messages, output_lines = zip(*results)
	for msg in messages:
		if msg is not None:
			raise ValueError(f"{msg}")

	with open(os.path.join(args.output_dir, f"{args.genome_id}.cogs.txt"), "wt") as cogs_out:
		print(
			*("cog", "query", "dbhit",	"bitscore", "identity",	"matches", "mismatches", "gaps", "query_start", "query_end", "dbhit_start",	"dbhit_end", "strand",	"specI_only:specI_cluster",	"combined_cf", "score_cf",),
			sep="\t", file=cogs_out, flush=True
		)
		for line in it.chain(*output_lines):
			print("\t".join(line), file=cogs_out)
			specis[line[14]] += 1 
		

		# for cog in COGS:
		# 	if args.marker_set == "motus" and not COGS[cog]:
		# 		continue
		# 	cog_file = os.path.join(cog_dir, f"{cog}.fna")
		# 	if os.path.isfile(cog_file):
		# 		with open(align_file, "wt") as aln_file, open(cog_file, "rt") as cog_in:
		# 			for line in cog_in:
		# 				if line[0] == ">":
		# 					line = re.sub(r"$", f"  # {cog} {args.genome_id}", line.strip())
		# 				print(line.strip(), file=aln_file)

		# 		_, cog_lines = call_mapseq(
		# 			align_file, args.cog_db, cog, speci_header=None,
		# 		)

		# 		if cog_lines is not None:
		# 			for line in cog_lines:
		# 				print("\t".join(line), file=cogs_out)
		# 				specis[line[14]] += 1
		# 			cogs_out.flush()


				# mapseq_pr = subprocess.Popen(
				# 	[
				# 		"mapseq",
				# 		align_file,
				# 		os.path.join(args.cog_db, f"{cog}.fna"),
				# 		os.path.join(args.cog_db, f"{cog}.specI.tax"),					
				# 	],
				# 	stdout=subprocess.PIPE, stderr=subprocess.PIPE,
				# )

				# out, err = mapseq_pr.communicate()

				# out = out.decode().strip().split("\n")
				# if speci_header is None:				
				# 	speci_header = [line.strip().split("\t") for line in out if line[0] == "#"]
				# 	speci_header[-1].insert(0, "cog")

				# speci_cog_d[cog] = [line.strip().split("\t") for line in out if line[0] != "#"]
				# for line in speci_cog_d[cog]:
				# 	line.insert(0, cog)
				
				# cat mapseq/\${bin_id}/speci/* | sed '3,\${ /^#/d }' > ${sample_id}/\${bin_id}.speci.assignments

				# break

	speci_out = open(os.path.join(args.output_dir, f"{args.genome_id}.specI.txt"), "wt")
	speci_status_out = open(os.path.join(args.output_dir, f"{args.genome_id}.specI.status"), "wt")
	seqfile = os.path.join(args.output_dir, f"{args.genome_id}.specI.ffn.gz")
	open(seqfile, "wt").close()

	with speci_out, speci_status_out:
		# specis = Counter()
		# if speci_header:
		# 	print(*("\t".join(line) for line in speci_header), sep="\n", file=cogs_out)
		# for cog in COGS:
		# 	cog_lines = speci_cog_d.get(cog)
		# 	if cog_lines is not None:
		# 		for line in cog_lines:
		# 			print("\t".join(line), file=cogs_out)
		# 			specis[line[14]] += 1
		
		
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
					print("Getting sequences from cluster:", first[0], seqfile, "...", end="")
					n_seqs = get_sequences_from_cluster(dbstr, first[0], seqfile)
					print(n_seqs)

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
