#!/usr/bin/env nextflow

/* Outlining the STAR process of creating a reference genome index */
process STAR_REFERENCE_INDEX { 
    label 'process_high'
    publishDir "${outputDir}/star", mode: "copy"

    container 'community.wave.seqera.io/library/star:2.7.11b--822039d47adf19a7'

	input:
        val outputDir
		path genome_fasta_files
        path gtf_file
        val starReferenceIndexArgs

	output:
        path reference_index, emit: reference_index
	
    script:
    """
    mkdir reference_index
    STAR --runThreadN $task.cpus --runMode genomeGenerate --genomeDir reference_index --genomeFastaFiles $genome_fasta_files --sjdbGTFfile $gtf_file $starReferenceIndexArgs
    """
}

/* Outlining the HISAT2 process of creating a reference genome index */
process HISAT2_REFERENCE_INDEX {
    label 'process_high'
    publishDir "${outputDir}/hisat2/reference_index", mode: "copy"

    container 'community.wave.seqera.io/library/hisat2:2.2.1--df34d2bb25ac6de5'

	input:
        val outputDir
        val hisat2_index_prefix
        path genome_fasta_files
        val hisat2ReferenceIndexArgs

	output:
        path "${hisat2_index_prefix}*", emit: hisat2_index_files
	
    script:
    """
    hisat2-build $genome_fasta_files $hisat2_index_prefix -p $task.cpus $hisat2ReferenceIndexArgs
    """
}

/* Outlining the SALMON_QUASI_MAPPING process of creating a reference genome index */
process SALMON_REFERENCE_INDEX {
    label 'process_high'
    publishDir "${outputDir}/salmon", mode: "copy"

    container 'community.wave.seqera.io/library/salmon:1.10.3--fcd0755dd8abb423'

	input:
        val outputDir
        path transcripts_file
        val kmer_size

	output:
        path "reference_index", emit: reference_index
	
    script:
    """
    mkdir reference_index
    salmon index -t $transcripts_file -k $kmer_size -i reference_index
    """
}

/* Outlining the KALLISTO process of creating a reference genome index */
process KALLISTO_REFERENCE_INDEX {
    label 'process_high'
    publishDir "${outputDir}/kallisto/reference_index", mode: "copy"

    container 'community.wave.seqera.io/library/kallisto:0.51.1--b63691b6841c7a52'

	input:
        val outputDir
        path transcripts_file
        val kallistoIndexArgs

	output:
        path "kallisto_index.idx", emit: kallisto_reference_index
	
    script:
    """
    kallisto index -i kallisto_index.idx $transcripts_file $kallistoIndexArgs
    """
}
