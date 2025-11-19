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

    label 'majiq_build'
    publishDir "${outputDir}/majiq/build", mode: "copy"

    container 'mcfonsecalab/majiq'

	input:
        path config
        path gff_file
        val outputDir

	output:
        path "*.majiq", emit: samples_splice_graphs
        path "*.sj", emit: samples_splice_junctions
        path "majiq.log", emit: majiq_log_file
        path "splicegraph.sql", emit: splicegraph_file
    
    script:
    """
    majiq build -c $config $gff_file -o .
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

	output:
        val "${grouped_file_names[0]}", emit: sample_group
        path "${grouped_file_names[0]}.psi.tsv", emit: majiq_tsv_file
        path "${grouped_file_names[0]}.psi.voila", emit: majiq_voila_file
        path "psi_majiq.log", emit: majiq_log_file
    
    script:
    majiq_file_names = grouped_file_names[1].collect { it + '.majiq' }.join(' ')
    """
    majiq psi $majiq_file_names -o . -n ${grouped_file_names[0]} -j $task.cpus
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

	output:
        val "${grouped_files_pairs[0]}_v_${grouped_files_pairs[2]}", emit: paired_samples_name
        path "${grouped_files_pairs[0]}-${grouped_files_pairs[2]}.deltapsi.tsv", emit: majiq_tsv_file
        path "${grouped_files_pairs[0]}-${grouped_files_pairs[2]}.deltapsi.voila", emit: majiq_voila_file
        path "deltapsi_majiq.log", emit: majiq_log_file
    
    script:
    majiq_file_names_1 = grouped_files_pairs[1].collect { it + '.majiq' }.join(' ')
    majiq_file_names_2 = grouped_files_pairs[3].collect { it + '.majiq' }.join(' ')
    """
    majiq deltapsi -grp1 $majiq_file_names_1 -grp2 $majiq_file_names_2 -o . -n ${grouped_files_pairs[0]} ${grouped_files_pairs[2]} -j $task.cpus
    """
}

/* Outlining the VOILA psi process */
process VOILA_PSI {

    label 'voila_psi'
    tag "${sample_group}"
    publishDir "${outputDir}/voila/voila_psi/${sample_group}_psi", mode: "copy"

    container 'mcfonsecalab/majiq'

	input:
        val sample_group
        path majiq_psi_files_parent
        path splicegraph_file
        val outputDir

	output:
        path "*"
    
    script:
    """
    voila modulize -d . $splicegraph_file $majiq_psi_files_parent -j $task.cpus
    """
}

/* Outlining the VOILA delta psi process */
process VOILA_DELTA_PSI {

    label 'voila_delta_psi'
    tag "${paired_samples_name}"
    publishDir "${outputDir}/voila/voila_delta_psi/${paired_samples_name}", mode: "copy"

    container 'mcfonsecalab/majiq'

	input:
        val paired_samples_name
        path majiq_delta_psi_files_parent
        path splicegraph_file
        val outputDir

	output:
        path "*"
    
    script:
    """
    voila modulize -d . $splicegraph_file $majiq_delta_psi_files_parent -j $task.cpus
    """
}

/* Outlining the SUPPA2 generate event annotations process */
process SUPPA2_GENERATE_EVENT_ANNOTATIONS {
    memory '7.6 GB'
    cpus 2

    label 'suppa2_generate_event_annotations'
    publishDir "${outputDir}/suppa2/event_annotations", mode: "copy"

    container 'naotokubota/suppa:2.3'

	input:
        path gtf_file
        val outputDir

	output:
        path "*.gtf", emit: gtf_files
        path "*.ioe", emit: ioe_files
        path "all_events.ioe", emit: ioe_file
    
    script:
    """
    python /SUPPA-2.3/suppa.py generateEvents -i $gtf_file -p -f ioe -e SE SS MX RI FL -o IOE
    # Put all the ioe events in the same file
    awk '
        FNR==1 && NR!=1 { while (/^<header>/) getline; }
        1 {print}
    ' *.ioe > all_events.ioe
    """
}

/* Outlining the SUPPA2 calculate psi for events process */
process SUPPA2_CALCULATE_EVENTS_PSI {
    memory '7.6 GB'
    cpus 2

    label 'suppa2_calculate_events_psi'
    publishDir "${outputDir}/suppa2/suppa2_psi", mode: "copy"

    container 'suppa2_with_r'

	input:
        val sample_ids_list  // e.g. ["AIY_1", "AIY_2"]
        val quants_files_list  // e.g. ["/path/to/quant.sf", "/path/to/another/quant.sf"]
        path all_ioe_events_file
        val tpm_column
        val outputDir

	output:
        path "samples_tpm.txt", emit: tpm_file
        path "all_samples.psi", emit: psi_file
    
    script:
    quant_files = quants_files_list.join(' ')
    """
    # Extract the TPM values from the 4th column of the salmon output
    python /SUPPA-2.3/multipleFieldSelection.py -i $quant_files -k 1 -f $tpm_column -o samples_tpm.txt
    # Get the psi values
    python /SUPPA-2.3/suppa.py psiPerEvent -i $all_ioe_events_file -e samples_tpm.txt -o all_samples
    """
}

/* Outlining the SUPPA2 splitting files before doing differential splicing analysis process */
process SUPPA2_SPLIT_FILES {

    label 'suppa2_split_files'
    tag "${grouped_files[0]}"
    publishDir "${outputDir}/suppa2/suppa2_psi/${grouped_files[0]}", mode: "copy"

    container 'suppa2_with_r'

	input:
        path r_file
        val grouped_files  // e.g.: "AIY", [AIY_1_sorted, AIY_2_sorted]
        path tpm_file
        path psi_file
        val outputDir

	output:
        tuple val("${grouped_files[0]}"), path("${grouped_files[0]}_events.psi"), path("${grouped_files[0]}.tpm"), emit: individual_sample_files
    
    script:
    group_ids = grouped_files[1].join(',')
    """
    /usr/bin/Rscript $r_file $tpm_file $group_ids ./${grouped_files[0]}.tpm
    /usr/bin/Rscript $r_file $psi_file $group_ids ./${grouped_files[0]}_events.psi
    """
}

/* Outlining the SUPPA2 calculate delta psi for events process */
process SUPPA2_CALCULATE_EVENTS_DELTA_PSI {
    memory '7.6 GB'
    cpus 2

    label 'suppa2_calculate_events_delta_psi'
    tag "${paired_files[0]}_v_${paired_files[3]}"
    publishDir "${outputDir}/suppa2/suppa2_delta_psi/${paired_files[0]}_v_${paired_files[3]}", mode: "copy"

    container 'suppa2_with_r'

	input:
        val paired_files  // e.g. ["AIY", "AIY_events.psi", "AIY.tpm", "ASK", "ASK_events.psi", "ASK.tpm"]
        path all_ioe_events_file
        val outputDir

	output:
        path "${paired_files[0]}_v_${paired_files[3]}_diffsplice.dpsi"
        path "${paired_files[0]}_v_${paired_files[3]}_diffsplice.psivec"
    
    script:
    """
    python /SUPPA-2.3/suppa.py diffSplice -m empirical -gc -i $all_ioe_events_file -p ${paired_files[1]} ${paired_files[4]} -e ${paired_files[2]} ${paired_files[5]} -o ${paired_files[0]}_v_${paired_files[3]}_diffsplice
    """
}
