#!/usr/bin/env nextflow

/* Outlining bowtie2's creating contaminant index process */
process BOWTIE2_CONTAMINANT_INDEX {
    label 'process_medium'
    publishDir "${outputDir}/bowtie2/contaminants_index", mode: "copy"

    container 'biocontainers/bowtie2:v2.4.1_cv1'

	input:
		path contaminants_files  // contaminant fasta files in a list, e.g. ['ncrna.fasta', 'trna.fasta']
        val outputDir

	output:
        path "contaminants.fa", emit: combined_contaminants_file
        path "contaminants_index.*", emit: contaminants_index_files
	
    script:
        """
        zcat -f ${contaminants_files.join(' ')} > contaminants.fa
        
        bowtie2-build contaminants.fa contaminants_index
        """
}

/* Outlining bowtie2's remove contaminants process */
process BOWTIE2_REMOVE_CONTAMINANTS {
    label 'process_medium'
    tag "${sample_id}"
    publishDir "${outputDir}/bowtie2/filtered_samples", mode: "copy"

    container 'biocontainers/bowtie2:v2.4.1_cv1'

	input:
        tuple val(sample_id), val(sample_group), path(reads)
		path contaminants_index_files
        val outputDir

	output:
        // For both single and paired-end reads
        tuple val(sample_id), val(sample_group), path("${sample_id}_clean*.fastq", arity: '1..*'), emit: cleaned_samples
	
    script:
    if (params.paired_end) {
        """
        bowtie2 -x contaminants_index -1 ${reads[0]} -2 ${reads[1]} --very-sensitive -p $task.cpus --un-conc ${sample_id}_clean.fastq -S /dev/null

        mv ${sample_id}_clean.1.fastq ${sample_id}_clean_1.fastq
        mv ${sample_id}_clean.2.fastq ${sample_id}_clean_2.fastq
        """
    } else {
        """
        bowtie2 -x contaminants_index -U ${reads[0]} --very-sensitive -p $task.cpus --un ${sample_id}_clean.fastq -S /dev/null
        """
    }
}
