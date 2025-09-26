#!/usr/bin/env nextflow

/* Outlining the SAM_TO_BAM process for converting sam output files from alignment to bam files */
process SAM_TO_BAM {
    label 'convert_sam_to_bam'
    tag "$sample_id"
    publishDir "${outputDir}/samtools/converted_to_bam_files", mode: "copy"

    container 'community.wave.seqera.io/library/samtools:1.22.1--eccb42ff8fb55509'

	input:
        val outputDir
        tuple val(sample_id), path(sam_file)

	output:
        val sample_id, emit: sample_id
        path "${sample_id}_Aligned.out.bam", emit: bam_file
	
    script:
    """
    samtools view -bS $sam_file > ${sample_id}_Aligned.out.bam
    """
}

/* Outlining the SORT_AND_INDEX_BAM process for sorting then indexing bam files before counting reads */
process SORT_AND_INDEX_BAM {
    label 'sort_and_index_bam_files'
    tag "$sample_id"
    publishDir "${outputDir}/samtools/sorted_and_indexed_bam_files", mode: "copy"

    container 'community.wave.seqera.io/library/samtools:1.22.1--eccb42ff8fb55509'

	input:
        tuple val(sample_id), path(bam_file)
        val outputDir

	output:
        val sample_id, emit: sample_id
        path "${sample_id}_sorted.bam", emit: sorted_bam_file
        path "${sample_id}_sorted.bam.bai", emit: indexed_bam_file
	
    script:
    """
    samtools sort -@ $task.cpus $bam_file -o ${sample_id}_sorted.bam
    samtools index -@ $task.cpus ${sample_id}_sorted.bam
    """
}
