import os
import re
import subprocess


def call_mapseq(align_file, cog_db, cog, threads=4, speci_header=None):
	mapseq_pr = subprocess.run(
		[
			"mapseq",
			"-nthreads", f"{threads}",
			align_file,
			os.path.join(cog_db, f"{cog}.fna"),
			os.path.join(cog_db, f"{cog}.specI.tax"),							
		],
		stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
	)

	# out, err = mapseq_pr.communicate()

	msg, speci_header, speci_cog = None, None, None
	if mapseq_pr.returncode != 0:
		msg = mapseq_pr.stdout
	else:
		out = mapseq_pr.stdout.decode().strip().split("\n") #out.decode().strip().split("\n")

		if speci_header is None:				
			speci_header = [line.strip().split("\t") for line in out if line[0] == "#"]
			speci_header[-1].insert(0, "cog")

		speci_cog = [[cog] + (line.strip().split("\t")) for line in out if line[0] != "#"]	

	return msg, speci_header, speci_cog


def task(cog_file, cog, genome_id, cog_db, threads):
	with open(cog_file + ".align", "wt") as aln_file, open(cog_file, "rt") as cog_in:
		for line in cog_in:
			if line[0] == ">":
				line = re.sub(r"$", f"  # {cog} {genome_id}", line.strip())
			print(line.strip(), file=aln_file)
	
	msg, _, cog_lines = call_mapseq(cog_file + ".align", cog_db, cog, threads=threads)
	return msg, cog_lines