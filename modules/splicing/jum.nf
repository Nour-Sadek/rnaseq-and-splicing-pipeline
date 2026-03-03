#!/usr/bin/env nextflow

/* Outlining the JUM_A.sh process */
process JUM_A {
    memory '7.6 GB'
    cpus 2

    label 'jum_A'
    tag "${grouped_files_pairs[0]}_v_${grouped_files_pairs[2]}"
    publishDir "${outputDir}/jum/${grouped_files_pairs[0]}_v_${grouped_files_pairs[2]}", mode: "copy"

    container 'jum:3.0.0'

	input:
        val grouped_files_pairs  // e.g.: "AIY", [AIY_1_Aligned.out.sam, AIY_2_Aligned.out.sam], "ASK", [ASK_1_Aligned.out.sam, ASK_2_Aligned.out.sam]
        val outputDir
        path aligned_out_files
        path sj_out_tab_files
        path bam_sorted_files
        // val jumArgs

	output:
        val "${grouped_files_pairs[0]}", emit: condition_1_name
        val "${grouped_files_pairs[1]}", emit: condition_2_name
        path "JUM_diff" 
	
    script:
    condition_1_file_names = grouped_files_pairs[1].collect { it.name.replaceFirst(/Aligned\.out\.sam$/, '') }.join(',')
    condition_2_file_names = grouped_files_pairs[3].collect { it.name.replaceFirst(/Aligned\.out\.sam$/, '') }.join(',')
    """
    bash /opt/JUM/JUM_A.sh --Folder /opt/JUM --JuncThreshold 5 --Condition1_fileNum_threshold 2 --Condition2_fileNum_threshold 2 --IRthreshold 5 --Readlength 100 --Thread $task.cpus --Condition1SampleName $condition_1_file_names --Condition2SampleName $condition_2_file_names
    """
}

process R_SCRIPT_JUM {
    label 'r_script_jum'
    tag "${grouped_files_pairs[0]}_v_${grouped_files_pairs[2]}"
    publishDir "${outputDir}/jum/${grouped_files_pairs[0]}_v_${grouped_files_pairs[2]}", mode: "copy"

    container 'jum:3.0.0'

    input:
        val grouped_files_pairs
        val outputDir
    
    output:
        path "AS_differential.txt"
    
    script:
    def rows = []

    // step through pairs: condition, file list
    for (int i = 0; i < grouped_files_pairs.size(); i += 2) {
        def condition = grouped_files_pairs[i]
        def files     = grouped_files_pairs[i + 1]

        files.each { f ->
            def sample = f.getName().replaceFirst(/_Aligned\.out\.sam$/, '')
            rows << [sample, condition]
        }   
    }

    // sort alphabetically by sample name
    rows = rows.sort { it[0] }

    // write output file
    def outFile = "experiment_design.txt"
    outFile = new File(outFile)
    outFile.text = "condition\n" + rows.collect { "${it[0]}\t${it[1]}" }.join("\n")
    """
    Rscript /opt/JUM/R_script_JUM.R experiment_design.txt > outputFile.Rout 2> errorFile.Rout
    """
}

/*
process JUM_B {
    memory '7.6 GB'
    cpus 2

    label 'jum_B'
    tag "${grouped_files_pairs[0]}_v_${grouped_files_pairs[2]}"
    publishDir "${outputDir}/jum/${grouped_files_pairs[0]}_v_${grouped_files_pairs[2]}", mode: "copy"

    container 'jum:3.0.0'

}
*/