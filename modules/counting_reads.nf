#!/usr/bin/env nextflow

/* Outlining the HTSEQ_COUNT reads quantification process */
process HTSEQ_COUNT {
    label 'htseq_count'
    tag "$sample_id"
    publishDir "${outputDir}/htseq_count", mode: "copy"

    container 'community.wave.seqera.io/library/htseq:2.0.9--4a65a9021e1142a5'

	input:
        tuple val(sample_id), val(sample_group), path(sorted_bam_file)
        val outputDir
        path gtf_file
        val htseqCountArgs

	output:
        val sample_id, emit: sample_id
		path "${sample_id}_htseq.txt", emit: htseq_file
	
    script:
    """
    htseq-count --nprocesses=${task.cpus} --format bam --order pos $htseqCountArgs $sorted_bam_file $gtf_file > ${sample_id}_htseq.txt
    """
}

/* Outlining the FEATURE_COUNTS reads quantification process */
process FEATURE_COUNTS {
    label 'feature_counts'
    tag "$sample_id"
    publishDir "${outputDir}/feature_counts", mode: "copy"

    container 'community.wave.seqera.io/library/subread:2.1.1--0ac4d7e46cd0c5d7'

	input:
        tuple val(sample_id), val(sample_group), path(sorted_bam_file)
        val outputDir
        path gtf_file
        val featureCountsArgs

	output:
		path "${sample_id}_feature.txt", emit: feature_file
	
    script:
    if (params.paired_end) {
        """
        featureCounts -a $gtf_file -o ${sample_id}_feature.txt -p -T $task.cpus $featureCountsArgs $sorted_bam_file
        """
    } else {
        """
        featureCounts -a $gtf_file -o ${sample_id}_feature.txt -T $task.cpus $featureCountsArgs $sorted_bam_file
        """
    }
}

/* Outlining the SALMON_ALIGNMENT_MODE reads quantification process (requires previous alignment to transcripts rather than genome) */
process SALMON_ALIGNMENT_MODE {
    label 'salmon_alignment_mode_counts'
    tag "$sample_id"
    publishDir "${outputDir}/salmon/salmon_alignment_mode_counts/${sample_id}", mode: "copy"

    container 'community.wave.seqera.io/library/salmon:1.10.3--fcd0755dd8abb423'

	input:
        tuple val(sample_id), val(sample_group), path(sorted_bam_file)
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
    publishDir "${outputDir}/salmon/salmon_quasi_mapping_mode_counts", mode: "copy"

    container 'community.wave.seqera.io/library/salmon:1.10.3--fcd0755dd8abb423'

	input:
        tuple val(sample_id), val(sample_group), path(reads)
        val outputDir
        path reference_index

	output:
        tuple val(sample_id), val(sample_group), path("${sample_id}/quant.sf"), emit: quants_file
        path "${sample_id}/aux_info", emit: aux_info
        path "${sample_id}/libParams", emit: lib_params
        path "${sample_id}/logs", emit: logs
        path "${sample_id}/cmd_info.json", emit: cmd_info
        path "${sample_id}/lib_format_counts.json", emit: lib_format_counts
	
    script:
    if (params.paired_end) {
        """
        mkdir $sample_id
        salmon quant -i $reference_index -l A -1 ${reads[0]} -2 ${reads[1]} -p $task.cpus --validateMappings --gcBias -o $sample_id
        """
    } else {
        """
        mkdir $sample_id
        salmon quant -i $reference_index -l A -r ${reads[0]} -p $task.cpus --validateMappings --gcBias -o $sample_id
        """
    }
}

/* Outlining the KALLISTO reads quantification process */
process KALLISTO {
    label 'kallisto'
    tag "$sample_id"
    publishDir "${outputDir}/kallisto/counts", mode: "copy"

    container 'community.wave.seqera.io/library/kallisto:0.51.1--b63691b6841c7a52'

	input:
        tuple val(sample_id), val(sample_group), path(reads)
        val outputDir
        path reference_index
        val num_bootstrap_samples

	output:
        tuple val(sample_id), val(sample_group), path("${sample_id}/abundance.tsv"), emit: quants_file
        path "${sample_id}/abundance.h5", emit: h5_file
        path "${sample_id}/run_info.json", emit: run_info
	
    script:
    if (params.paired_end) {
        """
        kallisto quant -i $reference_index -o ${sample_id} -t $task.cpus -b $num_bootstrap_samples ${reads[0]} ${reads[1]}
        """
    } else {
        """
        kallisto quant -i $reference_index -o ${sample_id} -t $task.cpus -b $num_bootstrap_samples --single -l $params.mean_fragment_length -s $params.sd_fragment_length ${reads[0]}
        """
    }
}

/* Outlining the RSEM alignment and reads quantification process */
process RSEM {
    memory '12 GB'
    cpus 2

    label 'rsem'
    tag "$sample_id"
    publishDir "${outputDir}/rsem/counts/${sample_id}", mode: "copy"

    container 'community.wave.seqera.io/library/rsem_star:f47af67e18e2d94b'

	input:
        tuple val(sample_id), val(sample_group), path(reads)
        val outputDir
        val rsem_index_prefix
        path rsem_index_files

	output:
		path "${sample_id}.stat", emit: stats_file
        path "${sample_id}.genes.results", emit: genes_results
        path "${sample_id}.isoforms.results", emit: isoforms_results
        path "${sample_id}.log", emit: log_file
        path "${sample_id}.transcript.bam", emit: transcript_bam_file
	
    script:
    if (params.paired_end) {
        """
        rsem-calculate-expression --paired-end --star -p $task.cpus ${reads[0]} ${reads[1]} reference_index/$rsem_index_prefix $sample_id
        """
    } else {
        """
        rsem-calculate-expression --star --strandedness $params.strandedness --fragment-length-mean $params.fragment_length_mean --fragment-length-sd $params.fragment_length_sd -p $task.cpus ${reads[0]} reference_index/$rsem_index_prefix $sample_id
        """
    }
}
