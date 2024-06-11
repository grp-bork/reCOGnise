params.recognise = [:]
params.recognise.marker_set = "motus"


process recognise {
	container "oras://registry.git.embl.de/schudoma/recognise-singularity/recognise-singularity:8b158eab"
	label "recognise"

	input:
	tuple val(speci), val(genome_id), path(genes), path(proteins)
	path(marker_genes_db)
	path(genes_db)
	path(db_credentials)
	
	output:
	tuple val(speci), val(genome_id), path("recognise/${genome_id}/${genome_id}.cogs.txt"), emit: cog_table
	tuple val(speci), val(genome_id), path("recognise/${genome_id}/${genome_id}.specI.txt"), emit: genome_speci
	tuple val(speci), val(genome_id), path("recognise/${genome_id}/${genome_id}.specI.status"), emit: speci_status
	tuple val(speci), val(genome_id), path("recognise/${genome_id}/${genome_id}.specI.ffn.gz"), emit: speci_sequences, optional: true
	tuple val(speci), val(genome_id), path("recognise/${genome_id}/${genome_id}.specI.status.OK"), emit: speci_status_ok, optional: true 

	script:
	"""
	recognise --marker_set ${params.recognise.marker_set} --genes ${genes} --proteins ${proteins} --cpus ${task.cpus} -o recognise/${genome_id} ${genome_id} \$(readlink ${marker_genes_db})
	"""

}