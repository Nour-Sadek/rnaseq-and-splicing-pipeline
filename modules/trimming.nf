#!/usr/bin/env nextflow

/* Outlining the Trimmomatic trimming process */
process TRIMMOMATIC {
    label 'trimmomatic'
    tag "$sample_id"
    publishDir "${outputDir}/trimmomatic/${sample_id}", mode: "copy"

    container 'community.wave.seqera.io/library/trimmomatic:0.40--0c25090769939729'

	input:
		tuple val(sample_id), path(read_1), path(read_2)
        path adapters_file
        val outputDir
        val trimmomaticArgs

	output:
    tuple val(sample_id), path("${sample_id}_fwd_trimmed.fastq"), path("${sample_id}_rev_trimmed.fastq"), emit: trimmed_samples
    path "${sample_id}_fwd_unpaired.fastq", emit: fwd_unpaired
    path "${sample_id}_rev_unpaired.fastq", emit: rev_unpaired
	
    script:
    """
    trimmomatic PE $read_1 $read_2 ${sample_id}_fwd_trimmed.fastq ${sample_id}_fwd_unpaired.fastq ${sample_id}_rev_trimmed.fastq ${sample_id}_rev_unpaired.fastq -threads $task.cpus ILLUMINACLIP:${adapters_file}:${trimmomaticArgs}
    """
}

/* Outlining the bbduk trimming process */
process BBDUK {
    label 'bbduk'
    tag "$sample_id"
    publishDir "${outputDir}/bbduk/${sample_id}", mode: "copy"

    container 'community.wave.seqera.io/library/bbmap:39.33--60639c9e1473b7a8'

	input:
		tuple val(sample_id), path(read_1), path(read_2)
        val outputDir
        val bbdukArgs

	output:
    tuple val(sample_id), path("${sample_id}_fwd_trimmed.fastq"), path("${sample_id}_rev_trimmed.fastq"), emit: trimmed_samples
    path("${sample_id}_unpaired.fastq"), emit: unpaired
	
    script:
    """
    bbduk.sh threads=$task.cpus ${bbdukArgs} in1=$read_1 in2=$read_2 out1="${sample_id}_fwd_trimmed.fastq" out2="${sample_id}_rev_trimmed.fastq" outs="${sample_id}_unpaired.fastq"
    """
}

/* Outlining the trim_galore trimming process */
process TRIM_GALORE {
    label 'trim_galore'
    tag "$sample_id"
    publishDir "${outputDir}/trim_galore/${sample_id}", mode: "copy"

    container 'community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18'

	input:
		tuple val(sample_id), path(read_1), path(read_2)
        val outputDir
        val trimGaloreArgs

	output:
    tuple val(sample_id), path("${read_1.simpleName}_val_1.fq"), path("${read_2.simpleName}_val_2.fq"), emit: trimmed_samples
    tuple path("${read_1.simpleName}_unpaired_1.fq"), path("${read_2.simpleName}_unpaired_2.fq"), emit: unpaired_reads
    tuple path("${read_1.simpleName}.fastq_trimming_report.txt"), path("${read_2.simpleName}.fastq_trimming_report.txt"), emit: trimming_reports
	
    script:
    """
    trim_galore --paired $read_1 $read_2 --cores $task.cpus --retain_unpaired ${trimGaloreArgs}
    """
}
