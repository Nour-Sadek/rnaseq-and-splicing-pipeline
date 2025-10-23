#!/usr/bin/env nextflow

/* Outlining the rMATS splciing process */
process rMATS {
    label 'rMats'
    tag "${grouped_bams[0][0]}_vs_${grouped_bams[1][0]}"
    publishDir "${outputDir}/rMats", mode: "copy"

    container 'mcfonsecalab/rmats'

	input:
        val grouped_bams
        path gtf_file
        val outputDir

	output:
        path "rmats_${grouped_bams[0][0]}_vs_${grouped_bams[1][0]}"
        path "rmats_temp_${grouped_bams[0][0]}_vs_${grouped_bams[1][0]}"

    script:
    """
    echo "${grouped_bams[0][1].join(',')}" > b1.txt
    echo "${grouped_bams[1][1].join(',')}" > b2.txt
    python /rmats-turbo/rmats.py --b1 b1.txt --b2 b2.txt --gtf $gtf_file --od rmats_${grouped_bams[0][0]}_vs_${grouped_bams[1][0]} --tmp rmats_temp_${grouped_bams[0][0]}_vs_${grouped_bams[1][0]} --readLength 40 --nthread $task.cpus
    
    """
}
