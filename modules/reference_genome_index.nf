#!/usr/bin/env nextflow

/* Outlining the STAR process of creating a reference genome index */
process STAR_REFERENCE_INDEX {
    label 'star_reference_index'
    publishDir "${outputDir}/star", mode: "copy"

    container 'community.wave.seqera.io/library/star:2.7.11b--822039d47adf19a7'

	input:
        val outputDir
		path genome_fasta_files
        path gtf_file
        val overhang
        val genomeSAindexNbases

	output:
        path reference_index, emit: reference_index
	
    script:
    """
    mkdir reference_index
    STAR --runThreadN $task.cpus --runMode genomeGenerate --genomeDir reference_index --genomeFastaFiles $genome_fasta_files --sjdbGTFfile $gtf_file --sjdbOverhang $overhang --genomeSAindexNbases $genomeSAindexNbases
    """
}

/* Outlining the HISAT2 process of creating a reference genome index */
process HISAT2_REFERENCE_INDEX {
    label 'hisat2_reference_index'
    publishDir "${outputDir}/hisat2/reference_index", mode: "copy"

    container 'community.wave.seqera.io/library/hisat2:2.2.1--df34d2bb25ac6de5'

	input:
        val outputDir
        val hisat2_index_prefix
		path genome_fasta_files

	output:
        path "${hisat2_index_prefix}*", emit: hisat2_index_files
        val "$hisat2_index_prefix", emit: hisat2_prefix_index
	
    script:
    """
    hisat2-build $genome_fasta_files $hisat2_index_prefix
    """
}

/* Outlining the HISAT2 process for improving the splicing alignment later */
process HISAT2_IMPROVE_SPLICE_ALIGNMENT {
    label 'hisat2_improve_splice_alignment'
    publishDir "${outputDir}/hisat2/improve_splice_alignment", mode: "copy"

    container 'community.wave.seqera.io/library/hisat2:2.2.1--df34d2bb25ac6de5'

	input:
        val outputDir
		path gtf_file

	output:
        path "genome.ss", emit: splice_sites
        path "genome.exon", emit: exons
	
    script:
    """
    hisat2_extract_splice_sites.py $gtf_file > genome.ss
    hisat2_extract_exons.py $gtf_file > genome.exon
    """
}

/* Outlining the SALMON_QUASI_MAPPING process of creating a reference genome index */
process SALMON_REFERENCE_INDEX {
    label 'salmon_reference_index'
    publishDir "${outputDir}/salmon", mode: "copy"

    container 'community.wave.seqera.io/library/salmon:1.10.3--fcd0755dd8abb423'

	input:
        val outputDir
        path transcripts_file
        val kmer_size

	output:
        path "reference_index", emit: reference_index
	
    script:
    """
    mkdir reference_index
    salmon index -t $transcripts_file -k $kmer_size -i reference_index
    """
}

/* Outlining the KALLISTO process of creating a reference genome index */
process KALLISTO_REFERENCE_INDEX {
    label 'kallisto_reference_index'
    publishDir "${outputDir}/kallisto/reference_index", mode: "copy"

    container 'community.wave.seqera.io/library/kallisto:0.51.1--b63691b6841c7a52'

	input:
        val outputDir
        path transcripts_file

	output:
        path "kallisto_index.idx", emit: kallisto_reference_index
	
    script:
    """
    kallisto index -i kallisto_index.idx $transcripts_file
    """
}

/* Outlining the RSEM process of creating a reference genome index */
process RSEM_REFERENCE_INDEX {
    label 'rsem_reference_index'
    publishDir "${outputDir}/rsem/reference_index", mode: "copy"

    container 'community.wave.seqera.io/library/rsem_star:f47af67e18e2d94b'

	input:
        val outputDir
        val rsem_index_prefix
        path genome_fasta_files
        path gtf_file

	output:
        val outputDir, emit: outputDir
        path "${rsem_index_prefix}.*", emit: rsem_index_files
        val "$rsem_index_prefix", emit: rsem_index_prefix
	
    script:
    """
    rsem-prepare-reference -gtf $gtf_file --star $genome_fasta_files $rsem_index_prefix
    """
}

/* Outlining the MINIMAP2 process of creating a reference genome index */
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