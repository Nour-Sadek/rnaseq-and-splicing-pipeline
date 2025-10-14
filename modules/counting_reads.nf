#!/usr/bin/env nextflow

/* Outlining the HTSEQ_COUNT reads quantification process */
process HTSEQ_COUNT {
    label 'htseq_count'
    tag "$sample_id"
    publishDir "${outputDir}/htseq_count", mode: "copy"

    container 'community.wave.seqera.io/library/htseq:2.0.9--4a65a9021e1142a5'

	input:
        tuple val(sample_id), path(sorted_bam_file)
        val outputDir
        path gtf_file

	output:
        val sample_id, emit: sample_id
		path "${sample_id}_htseq.txt", emit: htseq_file
	
    script:
    """
    htseq-count -f bam -r pos -s no -t exon -i gene_id -m union $sorted_bam_file $gtf_file > ${sample_id}_htseq.txt
    """
}

/* Outlining the FEATURE_COUNTS reads quantification process */
process FEATURE_COUNTS {
    label 'feature_counts'
    tag "$sample_id"
    publishDir "${outputDir}/feature_counts", mode: "copy"

    container 'community.wave.seqera.io/library/subread:2.1.1--0ac4d7e46cd0c5d7'

	input:
        tuple val(sample_id), path(sorted_bam_file)
        val outputDir
        path gtf_file

	output:
		path "${sample_id}_feature.txt", emit: feature_file
	
    script:
    """
    featureCounts -a $gtf_file -o ${sample_id}_feature.txt -t exon -g gene_id -s 0 -p -T ${task.cpus} $sorted_bam_file
    """
}

/* Outlining the SALMON_ALIGNMENT_MODE reads quantification process (requires previous alignment to transcripts rather than genome) */
process SALMON_ALIGNMENT_MODE {
    label 'salmon_alignment_mode_counts'
    tag "$sample_id"
    publishDir "${outputDir}/salmon/salmon_alignment_mode_counts/${sample_id}", mode: "copy"

    container 'community.wave.seqera.io/library/salmon:1.10.3--fcd0755dd8abb423'

	input:
        tuple val(sample_id), path(sorted_bam_file)
        val outputDir
        path transcripts_file

	output:
		path "*"
	
    script:
    """
    salmon quant -t $transcripts_file -l A -a $sorted_bam_file -p $task.cpus -o .
    """
}

/* Outlining the SALMON_QUASI_MAPPING_MODE reads quantification process */
process SALMON_QUASI_MAPPING_MODE {
    label 'salmon_quasi_mapping_mode_counts'
    tag "$sample_id"
    publishDir "${outputDir}/salmon/salmon_quasi_mapping_mode_counts/${sample_id}", mode: "copy"

    container 'community.wave.seqera.io/library/salmon:1.10.3--fcd0755dd8abb423'

	input:
        tuple val(sample_id), path(read_1), path(read_2)
        val outputDir
        path reference_index

	output:
		path "*"
	
    script:
    """
    salmon quant -i $reference_index -l A -1 $read_1 -2 $read_2 -p $task.cpus --validateMappings --gcBias -o .
    """
}

/* Outlining the KALLISTO reads quantification process */
process KALLISTO {
    label 'kallisto'
    tag "$sample_id"
    publishDir "${outputDir}/kallisto/counts/${sample_id}", mode: "copy"

    container 'community.wave.seqera.io/library/kallisto:0.51.1--b63691b6841c7a52'

	input:
        tuple val(sample_id), path(read_1), path(read_2)
        val outputDir
        path reference_index
        val num_bootstrap_samples

	output:
        path "abundance.h5"
        path "abundance.tsv"
        path "run_info.json"
	
    script:
    """
    kallisto quant -i $reference_index -o . -t $task.cpus -b $num_bootstrap_samples $read_1 $read_2
    """
}

/* Outlining the RSEM alignment and reads quantification process */
process RSEM {
    label 'rsem'
    tag "$sample_id"
    publishDir "${outputDir}/rsem/counts/${sample_id}", mode: "copy"

    container 'community.wave.seqera.io/library/rsem_star:f47af67e18e2d94b'

	input:
        tuple val(sample_id), path(read_1), path(read_2)
        val outputDir
        val rsem_index_prefix
        path rsem_index_files

	output:
		path "${sample_id}.stat"
        path "${sample_id}.genes.results"
        path "${sample_id}.isoforms.results"
        path "${sample_id}.log"
        path "${sample_id}.transcript.bam"
	
    script:
    """
    rsem-calculate-expression --paired-end --star -p $task.cpus $read_1 $read_2 $rsem_index_prefix $sample_id
    """
}
