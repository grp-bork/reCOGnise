params.recognise = [:]
params.recognise.marker_set = "motus"


process recognise_genome {
	container "ghcr.io/grp-bork/recognise:main"
	label "recognise"

	input:
	tuple val(speci), val(genome_id), path(genome)
	path(marker_genes_db)
	
	
	output:
	tuple val(genome_id), path("recognise/${genome_id}/${genome_id}.specI.status.OK"), emit: speci_status_ok, optional: true 
	tuple val(genome_id), path("recognise/${genome_id}/${genome_id}.cogs.txt"), emit: cog_table
	tuple val(genome_id), path("recognise/${genome_id}/${genome_id}.specI.txt"), emit: genome_speci
	tuple val(genome_id), path("recognise/${genome_id}/${genome_id}.specI.status"), emit: speci_status
	tuple val(genome_id), path("recognise/${genome_id}/${genome_id}.faa.gz"), emit: proteins
	tuple val(genome_id), path("recognise/${genome_id}/${genome_id}.ffn.gz"), emit: genes
	tuple val(genome_id), path("recognise/${genome_id}/${genome_id}.gff.gz"), emit: gff
	
	script:
	"""
	if [[ "${genome}" == *".gz" ]]; then
		gzip -dc ${genome} > genome_file
	else 
        cp -v ${genome} genome_file
    fi

	recognise --marker_set ${params.recognise.marker_set} --genome genome_file --cpus ${task.cpus} --with_gff -o recognise/${genome_id} ${genome_id} \$(readlink ${marker_genes_db})
	gzip -v recognise/${genome_id}/*.{faa,ffn,gff}

	rm -fv genome_file
	"""
}
