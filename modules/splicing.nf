#!/usr/bin/env nextflow

/* Outlining the rMATS individual splciing process */
process rMATS_INDIVIDUAL {
    memory '7.6 GB'
    cpus 2

    label 'rMats_individual'
    tag "${sample_bams[0]}"
    publishDir "${outputDir}/rMats/individual_analysis", mode: "copy"

    container 'community.wave.seqera.io/library/rmats:4.3.0--177f3a2035a879e5'

	input:
        val sample_bams
        val read_length
        path gtf_file
        val outputDir

	output:
        path "rmats_${sample_bams[0]}"
        path "rmats_temp_${sample_bams[0]}"
    
    script:
    """
    echo "${sample_bams[1].join(',')}" > b1.txt
    python /opt/conda/rMATS/rmats.py --b1 b1.txt --gtf $gtf_file --od rmats_${sample_bams[0]} --tmp rmats_temp_${sample_bams[0]} --readLength $read_length --nthread $task.cpus
    """
}

/* Outlining the rMATS differential splciing process */
process rMATS_DIFFERENTIAL {
    memory '7.6 GB'
    cpus 2
    
    label 'rMats_differential'
    tag "${grouped_bams[0][0]}_vs_${grouped_bams[1][0]}"
    publishDir "${outputDir}/rMats/differential_analysis", mode: "copy"

    container 'community.wave.seqera.io/library/rmats:4.3.0--177f3a2035a879e5'

	input:
        val grouped_bams
        val read_length
        path gtf_file
        val outputDir

	output:
        path "rmats_${grouped_bams[0][0]}_vs_${grouped_bams[1][0]}"
        path "rmats_temp_${grouped_bams[0][0]}_vs_${grouped_bams[1][0]}"

    script:
    """
    echo "${grouped_bams[0][1].join(',')}" > b1.txt
    echo "${grouped_bams[1][1].join(',')}" > b2.txt
    python /opt/conda/rMATS/rmats.py --b1 b1.txt --b2 b2.txt --gtf $gtf_file --od rmats_${grouped_bams[0][0]}_vs_${grouped_bams[1][0]} --tmp rmats_temp_${grouped_bams[0][0]}_vs_${grouped_bams[1][0]} --readLength $read_length --nthread $task.cpus
    """
}

/* Preparing the majiq config file */
process MAJIQ_CONFIG {
    label 'majiq_config'
    publishDir "${outputDir}/majiq/build", mode: "copy"

	input:
        val all_sample_bams
        val bam_dirs
        val outputDir

	output:
        path "config.ini", emit: config
    
    script:
    // Generate the [experiments] part of the config file
    experiments_lines = all_sample_bams.collect { g ->
        "${g[0]}=${g[1].join(',')}"
    }.join("\n")

    """
    echo "[info]" > config.ini
    echo "bamdirs=${bam_dirs.join(',')}" >> config.ini
    echo "genome=cel1\n" >> config.ini

    echo "[experiments]" >> config.ini
    echo "${experiments_lines}" >> config.ini
    """
}

/* Outlining the MAJIQ build process */
process MAJIQ_BUILD {
    memory '7.6 GB'
    cpus 2

    label 'majiq_build'
    publishDir "${outputDir}/majiq/build", mode: "copy"

    container 'mcfonsecalab/majiq'

	input:
        //tuple all_bam_files
        //tuple all_bam_bai_files
        path config
        path gff_file
        val outputDir

	output:
        path "*"
    
    script:
    """
    majiq build -c $config $gff_file -o .
    """
}


