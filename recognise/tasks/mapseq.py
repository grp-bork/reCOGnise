import os
import re
import subprocess


def call_mapseq(align_file, cog_db, cog, speci_header=None):
	mapseq_pr = subprocess.Popen(
		[
			"mapseq",
			align_file,
			os.path.join(cog_db, f"{cog}.fna"),
			os.path.join(cog_db, f"{cog}.specI.tax"),					
		],
		stdout=subprocess.PIPE, stderr=subprocess.PIPE,
	)

	out, err = mapseq_pr.communicate()

	msg, speci_header, speci_cog = None, None, None
	if mapseq_pr.returncode != 0:
		msg = err.decode().strip().split("\n")
	else:
		out = out.decode().strip().split("\n")

		if speci_header is None:				
			speci_header = [line.strip().split("\t") for line in out if line[0] == "#"]
			speci_header[-1].insert(0, "cog")

		speci_cog = [[cog] + (line.strip().split("\t")) for line in out if line[0] != "#"]	

	return msg, speci_header, speci_cog


def task(cog_file, cog, genome_id, cog_db):
	with open(cog_file + ".align", "wt") as aln_file, open(cog_file, "rt") as cog_in:
		for line in cog_in:
			if line[0] == ">":
				line = re.sub(r"$", f"  # {cog} {genome_id}", line.strip())
			print(line.strip(), file=aln_file)
	
	msg, _, cog_lines = call_mapseq(cog_file + ".align", cog_db, cog)
	return msg, cog_lines