#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { recognise_genome } from "./modules/recognise"


params.file_pattern = "**.{fna,fasta,fa,fna.gz,fasta.gz,fa.gz}"
print "PARAMS:\n" + params

def suffix_pattern = params.file_pattern.replaceAll(/\*/, "")


workflow {

    genomes_ch = Channel
        .fromPath(params.input_dir + "/" + params.file_pattern)
        .map { fasta ->
            def genome_id = fasta.name.replaceAll(/\.(fa|fasta|fna)(\.gz)?/, "")
            // def genome_id = fasta.name.split()
            return tuple(speci_tag, genome_id, fasta)
        }

    
    recognise_genome(genomes_ch, params.marker_db)
}

