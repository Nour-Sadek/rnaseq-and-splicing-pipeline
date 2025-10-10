#!/usr/bin/env nextflow

/* Outlining the GFFREAD transcriptome fasta generation utility process */
process GFFREAD {
    label 'gff_read'
    publishDir "${outputDir}/gffread", mode: "copy"

    container 'quay.io/biocontainers/gffread:0.12.7--h9a82719_0'

	input:
        val outputDir
        path genome_fasta_files
        path gtf_file

	output:
		path "transcripts.fa", emit: transcripts_fasta_file
	
    script:
    """
    gffread -w transcripts.fa -g $genome_fasta_files $gtf_file
    """
}