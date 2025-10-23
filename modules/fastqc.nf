#!/usr/bin/env nextflow

process FASTQC {
    
    label 'fastqc'
	tag "$sample_id"
    publishDir "${outputDir}", mode: "copy"

	container 'community.wave.seqera.io/library/fastqc:0.12.1--af7a5314d5015c29'

	input:
		tuple val(sample_id), val(sample_group), path(read_1), path(read_2)
        val(outputDir)

	output:
		path "*_fastqc.html", emit: html
		path "*_fastqc.zip", emit: zip
	
    script:
    """
    fastqc $read_1 $read_2 --threads ${task.cpus}
    """

}

