import argparse
import itertools as it
import logging
import multiprocessing as mp
import os
import pathlib
import sys

from collections import Counter
from contextlib import nullcontext

from .tasks.fetchmgs import call_fetch_mgs as fetchmgs
from .tasks.mapseq import task as mapseq
from .tasks.prodigal import call_prodigal as prodigal


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


def main():
    
    ap = argparse.ArgumentParser()
    ap.add_argument("genome_id", type=str)
    ap.add_argument("cog_db", type=str)
    ap.add_argument("--genes", type=str)
    ap.add_argument("--proteins", type=str)
    ap.add_argument("--genome", type=str)
    ap.add_argument("--cpus", type=int, default=4)
    ap.add_argument("--output_dir", "-o", type=str, default="recognise_out")
    ap.add_argument("--marker_set", type=str, choices=("full", "motus", "test"), default="motus")
    ap.add_argument("--with_gff", action="store_true")  # deprecated
    ap.add_argument("--min_markers", type=int, default=3)
    ap.add_argument("--min_clusters", type=int, default=2)
    ap.add_argument("--cluster_sizes", type=str)
    ap.add_argument("--create_workflow_sentinels", action="store_true")
    
    args = ap.parse_args()

    genome_present, genes_present, proteins_present = (
        f is not None and pathlib.Path(f).is_file()
        for f in (args.genome, args.genes, args.proteins)
    )
    genes, proteins = None, None

    accepted_clusters = None
    if args.cluster_sizes is not None and pathlib.Path(args.cluster_sizes).is_file():
        accepted_clusters = {
            line.strip().split("\t")[0]
            for line in open(args.cluster_sizes, "rt", encoding='utf-8')
            if int(line.strip().split("\t")[1]) >= args.min_clusters
        }
    
    output_dir = pathlib.Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    if genome_present:
        logger.info("Running prodigal...")
        if genes_present or proteins_present:
            raise ValueError("Please specify either a genome or a gene/protein set combination.")
        
        proteins = output_dir / f"{args.genome_id}.faa"
        genes = output_dir / f"{args.genome_id}.ffn"
        gff = output_dir / f"{args.genome_id}.gff"

        prodigal(args.genome, proteins, genes, gff)
        logger.info("prodigal finished.")

    elif genes_present and proteins_present:
        genes, proteins = args.genes, args.proteins
    elif genes_present:
        raise ValueError("Missing protein set, please specify with --proteins.")
    elif proteins_present:
        raise ValueError("Missing gene set, please specify with --genes.")

    cog_dir = output_dir / "cogs"
    
    logger.info("Running fetchMGs...")
    fetchmgs(proteins, genes, cog_dir, args.cpus)
    logger.info("fetchMGs finished.")

    specis = Counter()

    tasks = []
    for cog, is_motus_cog in COGS.items():
        if args.marker_set == "motus" and not is_motus_cog:
            continue
        if args.marker_set == "test" and len(tasks) == 3:
            break
        cog_file = cog_dir / f"{cog}.fna"
        if cog_file.is_file() and os.stat(cog_file).st_size:
            tasks.append((cog_file, cog, args.genome_id, pathlib.Path(args.cog_db), min(args.cpus, 4)))

    logger.info(
        f"Running {args.cpus // min(args.cpus, 4)} MAPseq processes on {len(tasks)} marker genes. "
        f"marker_set={args.marker_set}..."
    )

    with mp.Pool(args.cpus // min(args.cpus, 4)) as pool:
        # results = list(pool.starmap_async(task, tasks).get())
        results = [pool.apply_async(mapseq, (task,)) for task in tasks]
        # list(pool.apply_async(task, tasks).get())
        results = [res.get() for res in results]

    logger.info("MAPseq finished.")

    print(results)

    try:
        messages, output_lines = zip(*results)
    except ValueError:
        messages, output_lines = [], []

    for msg in messages:
        if msg is not None:
            raise ValueError(f"{msg}")

    with open(output_dir / f"{args.genome_id}.cogs.txt", "wt") as cogs_out:
        print(
            *("cog", "query", "dbhit",	"bitscore", "identity",	"matches", "mismatches", "gaps", "query_start", "query_end", "dbhit_start",	"dbhit_end", "strand",	"specI_only:specI_cluster",	"combined_cf", "score_cf",),
            sep="\t", file=cogs_out, flush=True
        )
        for line in it.chain(*output_lines):
            print("\t".join(line), file=cogs_out)
            specis[line[14]] += 1

    if args.create_workflow_sentinels:
        speci_out = open(args.output_dir / f"{args.genome_id}.specI.txt", "wt", encoding='utf-8',)
        speci_status_out = open(args.output_dir / f"{args.genome_id}.specI.status", "wt", encoding='utf-8',)
    else:
        speci_out, speci_status_out = nullcontext(), nullcontext()

    with speci_out, speci_status_out:
        speci_counts = specis.most_common()
        print(speci_counts)
        if not speci_counts:
            if args.create_workflow_sentinels:
                print("NO_MARKERS", file=speci_status_out)

            logger.warning("Could not find any markers. Aborting.")
        else:

            (speci_0, counts_0), *remaining = speci_counts
            speci_1, counts_1 = None, 0
            if remaining:
                (speci_1, counts_1), *remaining = remaining
            
            if counts_1 < counts_0 and 3 <= counts_0:
                if accepted_clusters and speci_0 not in accepted_clusters:
                    if args.create_workflow_sentinels:
                        print("SPECI_SIZE_INSUFFICIENT", file=speci_status_out)
                    logger.warning("specI cluster is too small. Aborting.")
                else:

                    logger.info("Found specI: %s" % speci_0)
                    if args.create_workflow_sentinels:
                        print(speci_0, file=speci_out)
                        print("OK", file=speci_status_out)
                    pathlib.Path(speci_status_out.name + ".OK").touch()

            else:
                if args.create_workflow_sentinels:
                    print("NO_CONSENSUS", file=speci_status_out)
                warn_params = (speci_0, counts_0, speci_1, counts_1,)
                logger.warning(
                    "Cannot determine consensus specI. first=%s (%s) second=%s (%s). Aborting." % warn_params
                )

    # """
    # for cog in \${specicogs[@]}; do
    # 		if [[ -s \${bin_id}_cogs/\${cog}.fna ]]
    # 		then
    # 			sed -i '/>/ s/\$/ '" # \${cog} \${bin_id}"'/' \${bin_id}_cogs/\${cog}.fna
    # 			((i=++i%6)) || wait
    # 			mapseq \${bin_id}_cogs/\${cog}.fna \${cogdir}\${cog}.fna \${cogdir}\${cog}.specI.tax | sed "s/#query/#cog\tquery/" | sed "1,2 ! s/^/\${cog}\t&/" > mapseq/\${bin_id}/speci/\${cog} &
    # 		else
    # 			touch mapseq/\${bin_id}/speci/\${cog}
    # 		fi
    # 	done
    # """
        
    
    

if __name__ == "__main__":
    main()

    # """
    # mkdir ${sample_id}
    # mkdir mapseq

    # rsync -av ${params.cog_ref_dir} \$TMPDIR
    # cogdir=\$TMPDIR/

    # specicogs=(COG0012 COG0016 COG0018 COG0048 COG0049 COG0052 COG0080 COG0081
    # 		   COG0085 COG0087 COG0088 COG0090 COG0091 COG0092 COG0093 COG0094
    # 		   COG0096 COG0097 COG0098 COG0099 COG0100 COG0102 COG0103 COG0124
    # 		   COG0172 COG0184 COG0185 COG0186 COG0197 COG0200 COG0201 COG0202
    # 		   COG0215 COG0256 COG0495 COG0522 COG0525 COG0533 COG0541 COG0552)
    # motu_cogs=(COG0012 COG0016 COG0018 COG0172 COG0215
    # 		   COG0495 COG0525 COG0533 COG0541 COG0552)
    # i=0

    # echo -n "bin" > ${sample_id}/${sample_id}.cog_count
    # for cog in \${specicogs[@]}; do
    # 	echo -n "\t"\${cog} >> ${sample_id}/${sample_id}.cog_count
    # done

    # for bin in bins/*
    # do
    #   bin_id=\${bin:5:\${#bin}-8}
    #   if [[ -s genecalls/\${bin_id}.extracted.fna  ]]
    #   then
    # 	fetchMG.pl -o \${bin_id}_cogs -t 5 -m extraction -d genecalls/\${bin_id}.extracted.fna genecalls/\${bin_id}.extracted.faa

    # 	mkdir mapseq/\${bin_id}
    # 	mkdir mapseq/\${bin_id}/motu
    # 	mkdir mapseq/\${bin_id}/speci

    # 	# SPECI
    # 	for cog in \${specicogs[@]}; do
    # 		if [[ -s \${bin_id}_cogs/\${cog}.fna ]]
    # 		then
    # 			sed -i '/>/ s/\$/ '" # \${cog} \${bin_id}"'/' \${bin_id}_cogs/\${cog}.fna
    # 			((i=++i%6)) || wait
    # 			mapseq \${bin_id}_cogs/\${cog}.fna \${cogdir}\${cog}.fna \${cogdir}\${cog}.specI.tax | sed "s/#query/#cog\tquery/" | sed "1,2 ! s/^/\${cog}\t&/" > mapseq/\${bin_id}/speci/\${cog} &
    # 		else
    # 			touch mapseq/\${bin_id}/speci/\${cog}
    # 		fi
    # 	done
    # 	wait
    # 	cat mapseq/\${bin_id}/speci/* | sed '3,\${ /^#/d }' > ${sample_id}/\${bin_id}.speci.assignments
    # 	i=0
    # 	# MOTU
    # 	for cog in \${motu_cogs[@]}; do
    # 		if [[ -s \${bin_id}_cogs/\${cog}.fna ]]
    # 		then
    # 			((i=++i%6)) || wait
    # 			mapseq \${bin_id}_cogs/\${cog}.fna \${cogdir}\${cog}.mOTU.fna \${cogdir}\${cog}.mOTU.tax | sed "s/#query/#cog\tquery/" | sed "1,2 ! s/^/\${cog}\t&/" > mapseq/\${bin_id}/motu/\${cog} &
    # 		else
    # 			touch mapseq/\${bin_id}/motu/\${cog}
    # 		fi
    # 	done
    # 	wait
    # 	cat mapseq/\${bin_id}/motu/* | sed '3,\${ /^#/d }' > ${sample_id}/\${bin_id}.motu.assignments

    # 	# SUMMARISE AND CLEANUP
    # 	cat \${bin_id}_cogs/*.fna >> ${sample_id}/${sample_id}.marker_genes.fna
    # 	taxonomy_annotation.py -m ${sample_id}/\${bin_id}.motu.assignments -s ${sample_id}/\${bin_id}.speci.assignments -o \${bin_id}.taxonomy_annotation -p \${bin_id}

    # 	# GET COG COUNTS
    # 	echo -n "\n"\${bin_id} >> ${sample_id}/${sample_id}.cog_count
    # 	for cog in \${specicogs[@]}; do
    # 		echo -n "\t"`grep -c \${cog} <(cat \${bin_id}_cogs/*.fna)` >> ${sample_id}/${sample_id}.cog_count
    # 	done
    # 	echo >> ${sample_id}/${sample_id}.cog_count
    #   else
    #   echo "No genecalls for \${bin_id}"
    #   fi
    # done

    # echo "name\ttaxon_level\ttaxonomic_assignment\tassignment_confidence\tconfident_count\tother_confident_count\tunconfident_count\tother_unconfident_count\tunique_marker_genes\ttotal_marker_genes\tchimerism_score_all\tchimerism_score_confident\tari_all\tari_confident" > ${sample_id}.taxonomy_annotation.tsv
    # cat *.taxonomy_annotation >> ${sample_id}.taxonomy_annotation.tsv

    # gzip ${sample_id}/${sample_id}.marker_genes.fna
    # mkdir ${sample_id}_assignments
    # mv ${sample_id}/*.assignments ${sample_id}_assignments/
    # tar -cvzf ${sample_id}/${sample_id}.assignments.tar.gz ${sample_id}_assignments/*.assignments
    # """
