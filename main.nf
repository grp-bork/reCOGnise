#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { species_recognition } from "./modules/recognise"


params.file_pattern = "**.fna"
print "PARAMS:\n" + params

def suffix_pattern = params.file_pattern.replaceAll(/\*/, "")


workflow {

    genomes_ch = Channel
        .fromPath(params.input_dir + "/" + params.file_pattern)
        .map { fasta ->
            def genome_id = fasta.name.replaceAll(suffix_pattern, "")
            return tuple(speci_tag, genome_id, fasta)
        }

    
    species_recognition(genomes_ch)
}

