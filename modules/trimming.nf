#!/usr/bin/env nextflow

/* Outlining the Trimmomatic trimming process */
process TRIMMOMATIC {
    label 'trimmomatic'
    tag "$sample_id"
    publishDir "${outputDir}/trimmomatic/${sample_id}", mode: "copy"

    container 'community.wave.seqera.io/library/trimmomatic:0.40--0c25090769939729'

	input:
		tuple val(sample_id), val(sample_group), path(reads)
        val outputDir
        val trimmomaticArgs

	output:
        // For both single and paired-end reads
        tuple val(sample_id), val(sample_group), path("${sample_id}*trimmed.fastq", arity: '1..*'), emit: trimmed_samples
        // For paired-end reads only
        path "${sample_id}_fwd_unpaired.fastq", emit: fwd_unpaired, optional: true
        path "${sample_id}_rev_unpaired.fastq", emit: rev_unpaired, optional: true
	
    script:
    if (params.paired_end) {
        """
        trimmomatic PE ${reads[0]} ${reads[1]} ${sample_id}_fwd_trimmed.fastq ${sample_id}_fwd_unpaired.fastq ${sample_id}_rev_trimmed.fastq ${sample_id}_rev_unpaired.fastq -threads $task.cpus ${trimmomaticArgs}
        """
    } else {
        """
        trimmomatic SE ${reads[0]} ${sample_id}_trimmed.fastq -threads $task.cpus ${trimmomaticArgs}
        """
    }
}

/* Outlining the bbduk trimming process */
process BBDUK {
    label 'bbduk'
    tag "$sample_id"
    publishDir "${outputDir}/bbduk/${sample_id}", mode: "copy"

    container 'community.wave.seqera.io/library/bbmap:39.33--60639c9e1473b7a8'

	input:
		tuple val(sample_id), val(sample_group), path(reads)
        val outputDir
        val bbdukArgs

	output:
        // For both single and paired-end reads
        tuple val(sample_id), val(sample_group), path("${sample_id}*trimmed.fastq", arity: '1..*'), emit: trimmed_samples
        // For paired-end reads only
        path("${sample_id}_unpaired.fastq"), emit: unpaired, optional: true
	
    script:
    if (params.paired_end) {
        """
        bbduk.sh -Xmx${task.memory.toGiga()}g threads=$task.cpus ${bbdukArgs} -da in1=${reads[0]} in2=${reads[1]} out1="${sample_id}_fwd_trimmed.fastq" out2="${sample_id}_rev_trimmed.fastq" outs="${sample_id}_unpaired.fastq"
        """
    } else {
        """
        bbduk.sh -Xmx${task.memory.toGiga()}g threads=$task.cpus ${bbdukArgs} -da in=${reads[0]} out="${sample_id}_trimmed.fastq"
        """
    }
}

/* Outlining the trim_galore trimming process */
process TRIM_GALORE {
    label 'trim_galore'
    tag "$sample_id"
    publishDir "${outputDir}/trim_galore/${sample_id}", mode: "copy"

    container 'community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18'

	input:
		tuple val(sample_id), val(sample_group), path(reads)
        val outputDir
        val trimGaloreArgs

	output:
        tuple val(sample_id), val(sample_group), path("*.fq", arity: '1..*'), emit: trimmed_samples
        tuple path("${read_1}_unpaired_1.fq"), path("${read_2}_unpaired_2.fq"), emit: unpaired_reads, optional: true
        tuple path("*.fastq_trimming_report.txt"), emit: trimming_reports
	
    script:
    if (params.paired_end) {
        read_1 = reads[0].simpleName
        read_2 = reads[1].simpleName
        """
        trim_galore --paired ${reads[0]} ${reads[1]} --cores $task.cpus --retain_unpaired ${trimGaloreArgs}
        """
    } else {
        """
        trim_galore ${reads[0]} --cores $task.cpus ${trimGaloreArgs}
        """
    }
}
