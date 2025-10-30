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
        val read_length
        path gtf_file
        val outputDir

	output:
        path "rmats_${sample_bams[0]}"
        path "rmats_temp_${sample_bams[0]}"
    
    script:
    """
    echo "${sample_bams[1].join(',')}" > b1.txt
    python /rmats-turbo/rmats.py --b1 b1.txt --gtf $gtf_file --od rmats_${sample_bams[0]} --tmp rmats_temp_${sample_bams[0]} --readLength $read_length --nthread $task.cpus
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
        val read_length
        path gtf_file
        val outputDir

	output:
        path "rmats_${sample_group_1}_vs_${sample_group_2}"
        path "rmats_temp_${sample_group_1}_vs_${sample_group_2}"

    script:
    """
    echo "${bam_files_1.join(',')}" > b1.txt
    echo "${bam_files_2.join(',')}" > b2.txt
    python /rmats-turbo/rmats.py --b1 b1.txt --b2 b2.txt --gtf $gtf_file --od rmats_${sample_group_1}_vs_${sample_group_2} --tmp rmats_temp_${sample_group_1}_vs_${sample_group_2} --readLength $read_length --nthread $task.cpus
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

/* Outlining the MAJIQ psi splicing process */
process MAJIQ_PSI {
    memory '7.6 GB'
    cpus 2

    label 'majiq_psi'
    tag "${sample_group}"
    publishDir "${outputDir}/majiq/majiq_psi/${sample_group}_psi", mode: "copy"

    container 'mcfonsecalab/majiq'

	input:
        tuple val(sample_group), val(file_names)  // e.g.: "AIY", [AIY_1_sorted, AIY_2_sorted]
        path build  // this folder contains the majiq files for the samples
        val outputDir

	output:
        path "${sample_group}.psi.tsv"
        path "${sample_group}.psi.voila"
        path "psi_majiq.log"
    
    script:
    majiq_file_names = file_names.collect { it + '.majiq' }.join(' ')
    """
    majiq psi $majiq_file_names -o . -n ${sample_group} -j $task.cpus
    """
}

/* Outlining the MAJIQ delta-psi splicing process */
process MAJIQ_DELTA_PSI {
    memory '7.6 GB'
    cpus 2

    label 'majiq_delta_psi'
    tag "${sample_group_1}_v_${sample_group_2}"
    publishDir "${outputDir}/majiq/majiq_delta_psi/${sample_group_1}_v_${sample_group_2}", mode: "copy"

    container 'mcfonsecalab/majiq'

	input:
        tuple val(sample_group_1), val(file_names_1), val(sample_group_2), val(file_names_2)  // e.g.: "AIY", [AIY_1_sorted, AIY_2_sorted], "ASK", [ASK_1_sorted, ASK_2_sorted]
        path build  // this folder contains the majiq files for the samples
        val outputDir

	output:
        path "${sample_group_1}-${sample_group_2}.deltapsi.tsv"
        path "${sample_group_1}-${sample_group_2}.deltapsi.voila"
        path "deltapsi_majiq.log"
    
    script:
    majiq_file_names_1 = file_names_1.collect { it + '.majiq' }.join(' ')
    majiq_file_names_2 = file_names_2.collect { it + '.majiq' }.join(' ')
    """
    majiq deltapsi -grp1 $majiq_file_names_1 -grp2 $majiq_file_names_2 -o . -n $sample_group_1 $sample_group_2 -j $task.cpus
    """
}
