#!/usr/bin/env nextflow

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
