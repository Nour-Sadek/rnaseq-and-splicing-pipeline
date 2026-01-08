#!/usr/bin/env nextflow

process FASTQC {
    
    label 'fastqc'
	tag "$sample_id"
    publishDir "${outputDir}", mode: "copy"

	container 'community.wave.seqera.io/library/fastqc:0.12.1--af7a5314d5015c29'

	input:
		tuple val(sample_id), val(sample_group), path(reads)  // reads could be either [read_1] for single-end reads or [read_1, read_2] for paired-end reads
        val(outputDir)

	output:
		path "*_fastqc.html", emit: html
		path "*_fastqc.zip", emit: zip
	
    script:
    if (params.paired_end) {
        """
        fastqc ${reads[0]} ${reads[1]} --threads ${task.cpus}
        """
    } else {
        """
        fastqc ${reads[0]} --threads ${task.cpus}
        """
    }

}

