#!/usr/bin/env nextflow

/* Outlining the rMATS individual splciing process */
process rMATS_INDIVIDUAL {
    memory '7.6 GB'
    cpus 2

    label 'rMats_individual'
    tag "${sample_bams[0]}"
    publishDir "${outputDir}/rMats/individual_analysis", mode: "copy"

    container 'mcfonsecalab/rmats'

	input:
        val sample_bams
        path gtf_file
        val outputDir
        val rmatsIndividualArgs

	output:
        path "rmats_${sample_bams[0]}"
        path "rmats_temp_${sample_bams[0]}"
    
    script:
    """
    echo "${sample_bams[1].join(',')}" > b1.txt
    python /rmats-turbo/rmats.py --b1 b1.txt --gtf $gtf_file --od rmats_${sample_bams[0]} --tmp rmats_temp_${sample_bams[0]} --nthread $task.cpus $rmatsIndividualArgs
    """
}

/* Outlining the rMATS differential splciing process */
process rMATS_DIFFERENTIAL {
    memory '7.6 GB'
    cpus 2
    
    label 'rMats_differential'
    tag "${sample_group_1}_v_${sample_group_2}"
    publishDir "${outputDir}/rMats/differential_analysis", mode: "copy"

    container 'mcfonsecalab/rmats'

	input:
        tuple val(sample_group_1), val(bam_files_1), val(sample_group_2), val(bam_files_2)
        path gtf_file
        val outputDir
        val rmatsDifferentialArgs

	output:
        path "rmats_${sample_group_1}_vs_${sample_group_2}"
        path "rmats_temp_${sample_group_1}_vs_${sample_group_2}"

    script:
    """
    echo "${bam_files_1.join(',')}" > b1.txt
    echo "${bam_files_2.join(',')}" > b2.txt
    python /rmats-turbo/rmats.py --b1 b1.txt --b2 b2.txt --gtf $gtf_file --od rmats_${sample_group_1}_vs_${sample_group_2} --tmp rmats_temp_${sample_group_1}_vs_${sample_group_2} --nthread $task.cpus $rmatsDifferentialArgs
    """
}
