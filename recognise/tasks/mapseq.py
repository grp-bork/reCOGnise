""" Module to call MAPseq """

import os
import re
import subprocess


def call_mapseq(align_file, cog_db, cog, threads=4, speci_header=None):
    """ Calls MAPseq """

    msg, speci_header, speci_cog = None, None, None

    try:
        mapseq_pr = subprocess.run(
            [
                "mapseq",
                "-nthreads", f"{threads}",
                align_file,
                os.path.join(cog_db, f"{cog}.fna"),
                os.path.join(cog_db, f"{cog}.specI.tax"),
            ],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            check=True, capture_output=True,
        )
    except subprocess.CalledProcessError as e:
        msg = f"MAPseq: {e.returncode}:\n{e.output}"
    else:
        out = mapseq_pr.stdout.decode().strip().split("\n")
        if speci_header is None:
            speci_header = [line.strip().split("\t") for line in out if line[0] == "#"]
            speci_header[-1].insert(0, "cog")

        speci_cog = [[cog] + (line.strip().split("\t")) for line in out if line[0] != "#"]

    with open(f"{align_file}.done", "wt", encoding="utf-8") as _out:
        ...

    return msg, speci_header, speci_cog


# def task(cog_file, cog, genome_id, cog_db, threads):
def task(argv):
    """ Preprocesses input fasta and calls MAPseq """
    cog_file, cog, genome_id, cog_db, threads = argv

    aln_file = open(cog_file + ".align", "wt", encoding="UTF-8",)
    cog_in = open(cog_file, "rt", encoding="UTF-8",)

    with aln_file, cog_in:
        for line in cog_in:
            if line[0] == ">":
                line = re.sub(r"$", f"  # {cog} {genome_id}", line.strip())
            print(line.strip(), file=aln_file)

    msg, _, cog_lines = call_mapseq(cog_file + ".align", cog_db, cog, threads=threads)
    return msg, cog_lines
