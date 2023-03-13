import gzip

import pymongo


def get_sequences_from_cluster(mongo_db_str, cluster_id, seqfile):

	client = pymongo.MongoClient(mongo_db_str,)
	fr_db = client["progenomes"]

	n_genes = 0
	files = []
	with gzip.open(seqfile, "wt") as genes_out:
		for record in fr_db.samples.find({'fr13_cluster': cluster_id}):
			genes_file = f"{record['analysis_path']}/ref_genome_called_genes/{record['sample_id']}.genes.fa.gz"
			files.append(genes_file)
		for genes_file in files:
			with gzip.open(genes_file, "rt") as genes_in:
				genes_raw = genes_in.read()
				n_genes += genes_raw.count(">")
				print(genes_raw, file=genes_out, end="" if genes_raw[-1] == "\n" else "\n")

	return n_genes
