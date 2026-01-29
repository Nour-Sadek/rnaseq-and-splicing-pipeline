#!/usr/bin/env nextflow

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

    label 'majiq_build'
    publishDir "${outputDir}/majiq/build", mode: "copy"

    container 'mcfonsecalab/majiq'

	input:
        path config
        path gff_file
        val outputDir
        val majiqBuildArgs

	output:
        path "*.majiq", emit: samples_splice_graphs
        path "*.sj", emit: samples_splice_junctions
        path "majiq.log", emit: majiq_log_file
        path "splicegraph.sql", emit: splicegraph_file
    
    script:
    """
    majiq build -c $config $gff_file -o . $majiqBuildArgs
    """
}

/* Outlining the MAJIQ psi splicing process */
process MAJIQ_PSI {

    label 'majiq_psi'
    tag "${grouped_file_names[0]}"
    publishDir "${outputDir}/majiq/majiq_psi/${grouped_file_names[0]}_psi", mode: "copy"

    container 'mcfonsecalab/majiq'

	input:
        val grouped_file_names  // e.g.: "AIY", [AIY_1_sorted, AIY_2_sorted]
        path build  // this folder contains the majiq files for the samples
        val outputDir
        val majiqPsiArgs

	output:
        val "${grouped_file_names[0]}", emit: sample_group
        path "${grouped_file_names[0]}.psi.tsv", emit: majiq_tsv_file
        path "${grouped_file_names[0]}.psi.voila", emit: majiq_voila_file
        path "psi_majiq.log", emit: majiq_log_file
    
    script:
    majiq_file_names = grouped_file_names[1].collect { it + '.majiq' }.join(' ')
    """
    majiq psi $majiq_file_names -o . -n ${grouped_file_names[0]} -j $task.cpus $majiqPsiArgs
    """
}

/* Outlining the MAJIQ delta-psi splicing process */
process MAJIQ_DELTA_PSI {

    label 'majiq_delta_psi'
    tag "${grouped_files_pairs[0]}_v_${grouped_files_pairs[2]}"
    publishDir "${outputDir}/majiq/majiq_delta_psi/${grouped_files_pairs[0]}_v_${grouped_files_pairs[2]}", mode: "copy"

    container 'mcfonsecalab/majiq'

	input:
        val grouped_files_pairs  // e.g.: "AIY", [AIY_1_sorted, AIY_2_sorted], "ASK", [ASK_1_sorted, ASK_2_sorted]
        path samples_splice_graphs  // this folder contains the majiq files for the samples
        val outputDir
        val majiqDeltaPsiArgs

	output:
        val "${grouped_files_pairs[0]}_v_${grouped_files_pairs[2]}", emit: paired_samples_name
        path "${grouped_files_pairs[0]}-${grouped_files_pairs[2]}.deltapsi.tsv", emit: majiq_tsv_file
        path "${grouped_files_pairs[0]}-${grouped_files_pairs[2]}.deltapsi.voila", emit: majiq_voila_file
        path "deltapsi_majiq.log", emit: majiq_log_file
    
    script:
    majiq_file_names_1 = grouped_files_pairs[1].collect { it + '.majiq' }.join(' ')
    majiq_file_names_2 = grouped_files_pairs[3].collect { it + '.majiq' }.join(' ')
    """
    majiq deltapsi -grp1 $majiq_file_names_1 -grp2 $majiq_file_names_2 -o . -n ${grouped_files_pairs[0]} ${grouped_files_pairs[2]} -j $task.cpus $majiqDeltaPsiArgs
    """
}
