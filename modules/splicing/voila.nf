#!/usr/bin/env nextflow

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
        val voilaModulizePsiArgs

	output:
        path "*"
    
    script:
    """
    voila modulize -d . $splicegraph_file $majiq_psi_files_parent -j $task.cpus $voilaModulizePsiArgs
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
        val voilaModulizeDeltaPsiArgs

	output:
        path "*"
    
    script:
    """
    voila modulize -d . $splicegraph_file $majiq_delta_psi_files_parent -j $task.cpus $voilaModulizeDeltaPsiArgs
    """
}
