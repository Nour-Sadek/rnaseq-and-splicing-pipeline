#!/usr/bin/env nextflow

/* Outlining the MINIMAP2 process of creating a reference genome index 
- previously in reference_genome_index.nf */
process MINIMAP2_REFERENCE_INDEX {
    label 'minimap2_reference_index'
    publishDir "${outputDir}/minimap2", mode: "copy"

    container 'community.wave.seqera.io/library/minimap2:2.30--dde6b0c5fbc82ebd'

	input:
        val outputDir
		path genome_fasta_files

	output:
        path reference_index.mmi, emit: reference_index
	
    script:
    """
    minimap2 -d reference_index.mmi $genome_fasta_files
    """
}

/* Outlining the MINIMAP2 alignment process 
- previously in aligning.nf */
process MINIMAP2 {
    memory '7 GB'
    cpus 2

    label 'minimap2'
    tag "$sample_id"
    publishDir "${outputDir}/minimap2/${sample_id}", mode: "copy"

    container 'community.wave.seqera.io/library/minimap2:2.30--dde6b0c5fbc82ebd'

	input:
        tuple val(sample_id), val(sample_group), path(reads)
        path reference_genome_index
        val outputDir
        val preset

	output:
        tuple val(sample_id), val(sample_group), path("${sample_id}_Aligned.sam"), emit: alignment_output
	
    script:
    if (params.paired_end) {
        """
        minimap2 -t $task.cpus -ax $preset $reference_genome_index ${reads[0]} ${reads[1]} > ${sample_id}_Aligned.sam
        """
    } else {
        """
        minimap2 -t $task.cpus -ax $preset $reference_genome_index ${reads[0]} > ${sample_id}_Aligned.sam
        """
    }
}

/* Outlining the SALMON_ALIGNMENT_MODE reads quantification process (requires previous alignment to transcripts rather than genome)
- previously in counting_reads.nf */
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

/* Outlining the RSEM process of creating a reference genome index
- previously in reference_genome_index.nf */
process RSEM_REFERENCE_INDEX {
    memory '12 GB'
    cpus 2

    label 'rsem_reference_index'
    publishDir "${outputDir}/rsem/reference_index", mode: "copy"

    container 'community.wave.seqera.io/library/rsem_star:f47af67e18e2d94b'

	input:
        val outputDir
        val rsem_index_prefix
        path genome_fasta_files
        path gtf_file
        val overhang
        val genomeSAindexNbases

	output:
        path reference_index, emit: rsem_index_files
	
    script:
    """
    mkdir reference_index
    rsem-prepare-reference -gtf $gtf_file $genome_fasta_files reference_index/$rsem_index_prefix
    STAR --runThreadN $task.cpus --runMode genomeGenerate --genomeDir reference_index --genomeFastaFiles $genome_fasta_files --sjdbGTFfile $gtf_file --sjdbOverhang $overhang --genomeSAindexNbases $genomeSAindexNbases
    """
}

/* Outlining the RSEM alignment and reads quantification process
- previously in counting_reads.nf */
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

/* Outlining the GFFREAD transcriptome fasta generation utility process
- previously in gff_utilities.nf */
process GFFREAD {
    label 'gff_read'
    publishDir "${outputDir}/gffread", mode: "copy"

    container 'quay.io/biocontainers/gffread:0.12.7--h9a82719_0'

	input:
        val outputDir
        path genome_fasta_files
        path gtf_file

	output:
		path "transcripts.fa", emit: transcripts_fasta_file
	
    script:
    """
    gffread -w transcripts.fa -g $genome_fasta_files $gtf_file
    """
}
