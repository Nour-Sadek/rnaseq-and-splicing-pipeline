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
        tuple val(sample_id), path("${sample_id}_Aligned.out.bam"), emit: alignment_output
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
    memory '7 GB'
    cpus 2

    label 'hisat2'
    tag "$sample_id"
    publishDir "${outputDir}/hisat2/${sample_id}", mode: "copy"

    container 'community.wave.seqera.io/library/hisat2:2.2.1--df34d2bb25ac6de5'

	input:
        tuple val(sample_id), path(read_1), path(read_2)
        val outputDir
        val hisat2_prefix_index
        path hisat2_index_files
        path splice_sites
        path exons

	output:
        tuple val(sample_id), path("${sample_id}_Aligned.out.sam"), emit: alignment_output
	
    script:
    """
    hisat2 -q -x $hisat2_prefix_index -ss $splice_sites -exon $exons -1 $read_1 -2 $read_2 -S ${sample_id}_Aligned.out.sam
    """
}

/* Outlining the MINIMAP2 alignment process */
process MINIMAP2 {
    memory '7 GB'
    cpus 2

    label 'minimap2'
    tag "$sample_id"
    publishDir "${outputDir}/minimap2/${sample_id}", mode: "copy"

    container 'community.wave.seqera.io/library/minimap2:2.30--dde6b0c5fbc82ebd'

	input:
        tuple val(sample_id), path(read_1), path(read_2)
        path reference_genome_index
        val outputDir
        val preset

	output:
        tuple val(sample_id), path("${sample_id}_Aligned.sam"), emit: alignment_output
	
    script:
    """
    minimap2 -t $task.cpus -ax $preset $reference_genome_index $read_1 $read_2 > ${sample_id}_Aligned.sam
    """
}