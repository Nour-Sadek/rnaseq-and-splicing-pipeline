#!/usr/bin/env nextflow

/* Outlining the STAR alignment process */
process STAR {
    label 'process_high'
    tag "$sample_id"
    publishDir "${outputDir}/star/${sample_id}", mode: "copy"

    container 'community.wave.seqera.io/library/star:2.7.11b--822039d47adf19a7'

	input:
        tuple val(sample_id), val(sample_group), path(reads)
        val outputDir
        path reference_genome_index
        val starArgs

	output:
        tuple val(sample_id), val(sample_group), path("${sample_id}_Aligned.out.sam"), emit: alignment_output
        path "${sample_id}_Aligned.out.sam", emit: sam_output
        path "${sample_id}_Log.final.out", emit: log_final
        path "${sample_id}_Log.out", emit: log
        path "${sample_id}_Log.progress.out", emit: log_progress
        path "${sample_id}_SJ.out.tab", emit: splicing_junctions
	
    script:
    if (params.paired_end) {
        """
        STAR --runThreadN $task.cpus --genomeDir $reference_genome_index --readFilesIn ${reads[0]} ${reads[1]} --outFileNamePrefix ${sample_id}_ $starArgs
        """
    } else {
        """
        STAR --runThreadN $task.cpus --genomeDir $reference_genome_index --readFilesIn ${reads[0]} --outFileNamePrefix ${sample_id}_ $starArgs
        """
    }
}

/* Outlining the HISAT2 alignment process */
process HISAT2 {
    label 'process_high'
    tag "$sample_id"
    publishDir "${outputDir}/hisat2/${sample_id}", mode: "copy"

    container 'community.wave.seqera.io/library/hisat2:2.2.1--df34d2bb25ac6de5'

	input:
        tuple val(sample_id), val(sample_group), path(reads)
        val outputDir
        val hisat2_prefix_index
        path hisat2_index_files
        val hisat2Args

	output:
        tuple val(sample_id), val(sample_group), path("${sample_id}_Aligned.out.sam"), emit: alignment_output
	
    script:
    if (params.paired_end) {
        """
        hisat2 -q -p $task.cpus -x $hisat2_prefix_index -1 ${reads[0]} -2 ${reads[1]} -S ${sample_id}_Aligned.out.sam $hisat2Args
        """
    } else {
        """
        hisat2 -q -p $task.cpus -x $hisat2_prefix_index -U ${reads[0]} -S ${sample_id}_Aligned.out.sam $hisat2Args
        """
    }
}

