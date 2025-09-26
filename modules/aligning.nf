#!/usr/bin/env nextflow

/* Outlining the STAR alignment process */
process STAR {
    memory '7 GB'
    cpus 2

    label 'star'
    tag "$sample_id"
    publishDir "${outputDir}/star/${sample_id}", mode: "copy"

    container 'community.wave.seqera.io/library/star:2.7.11b--822039d47adf19a7'

	input:
        tuple val(sample_id), path(read_1), path(read_2)
        val outputDir
        path reference_genome_index
        val filterMatch

	output:
        val sample_id, emit: sample_id
		path "${sample_id}_Aligned.out.bam", emit: bam_file
        path "${sample_id}_Log.final.out", emit: log_final
        path "${sample_id}_Log.out", emit: log
        path "${sample_id}_Log.progress.out", emit: log_progress
        path "${sample_id}_SJ.out.tab", emit: splicing_junctions
	
    script:
    """
    STAR --runThreadN $task.cpus --genomeDir $reference_genome_index --readFilesIn $read_1 $read_2 --outFileNamePrefix ${sample_id}_ --outSAMtype BAM Unsorted --outFilterMatchNminOverLread $filterMatch
    """
}

/* Outlining the HISAT2 alignment process */
process HISAT2 {
    label 'hisat2'
    tag "$sample_id"
    publishDir "${outputDir}/hisat2/${sample_id}", mode: "copy"

    container 'community.wave.seqera.io/library/hisat2:2.2.1--df34d2bb25ac6de5'

	input:
        tuple val(sample_id), path(read_1), path(read_2)
        val outputDir
        path reference_genome_index
        path splice_sites
        path exons

	output:
        val sample_id, emit: sample_id
		path "${sample_id}_Aligned.out.sam", emit: sam_file
	
    script:
    """
    hisat2 -x $reference_genome_index -ss $splice_sites -exon $exons -1 $read_1 -2 $read_2 -S ${sample_id}_Aligned.out.sam
    """
}
