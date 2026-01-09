#!/usr/bin/env nextflow

/* Outlining the WHIPPET creating index process */
process WHIPPET_INDEX {
    memory '7.6 GB'
    cpus 2

    label 'whippet_index'
    publishDir "${outputDir}/whippet/index", mode: "copy"

    container 'naotokubota/whippet:1.6.1'

	input:
        path genome_fasta_file
        path gtf_file
        val outputDir

	output:
        path "whippet_index.exons.tab.gz", emit: exons
        path "whippet_index.jls", emit: whippet_index
    
    script:
    """
    /usr/local/julia/bin/julia /Whippet.jl/bin/whippet-index.jl --fasta $genome_fasta_file --gtf $gtf_file --index whippet_index
    """
}

/* Outlining the WHIPPET quantifying FASTQ file process */
process WHIPPET_QUANT {
    memory '7.6 GB'
    cpus 2

    label 'whippet_quant'
    tag "${sample_id}"
    publishDir "${outputDir}/whippet/quant/${sample_id}", mode: "copy"

    container 'naotokubota/whippet:1.6.1'

	input:
        tuple val(sample_id), val(sample_group), path(reads)
        path whippet_index
        val outputDir

	output:
        tuple val(sample_id), val(sample_group), path("${sample_id}.psi.gz"), emit: sample_psi_file
        path "${sample_id}.gene.tpm.gz", emit: gene_tpm_file
        path "${sample_id}.isoform.tpm.gz", emit: isoform_tpm_file
        path "${sample_id}.jnc.gz", emit: junctions_file
        path "${sample_id}.map.gz", emit: map_file
    
    script:
    read_1 = reads[0]
    if (params.paired_end) {
        read_2 = reads[1]
        """
        # Fix the fatsq files by letting them follow the standard 4 line format
        awk 'NR % 4 == 3 { print "+"; next } { print }' "$read_1" > fixed_R1.fastq
        awk 'NR % 4 == 3 { print "+"; next } { print }' "$read_2" > fixed_R2.fastq
        /usr/local/julia/bin/julia /Whippet.jl/bin/whippet-quant.jl fixed_R1.fastq fixed_R2.fastq -x $whippet_index -o $sample_id
        """
    } else {
        """
        # Fix the fatsq files by letting them follow the standard 4 line format
        awk 'NR % 4 == 3 { print "+"; next } { print }' "$read_1" > fixed_R1.fastq
        /usr/local/julia/bin/julia /Whippet.jl/bin/whippet-quant.jl fixed_R1.fastq -x $whippet_index -o $sample_id
        """
    }
}

/* Outlining the WHIPPET delta psi process */
process WHIPPET_DELTA {
    memory '7.6 GB'
    cpus 2

    label 'whippet_delta'
    tag "${grouped_files_pairs[0]}_v_${grouped_files_pairs[2]}"
    publishDir "${outputDir}/whippet/delta_psi/${grouped_files_pairs[0]}_v_${grouped_files_pairs[2]}", mode: "copy"

    container 'naotokubota/whippet:1.6.1'

	input:
        val grouped_files_pairs  // e.g.: "AIY", [AIY_1.psi.gz, AIY_2.psi.gz], "ASK", [ASK_1.psi.gz, ASK_2.psi.gz]
        path psi_files
        val outputDir

	output:
        path "${grouped_files_pairs[0]}_v_${grouped_files_pairs[2]}.diff.gz"
    
    script:
    group_one = grouped_files_pairs[1].join(',')
    group_two = grouped_files_pairs[3].join(',')
    """
    /usr/local/julia/bin/julia /Whippet.jl/bin/whippet-delta.jl -a $group_one -b $group_two -o ${grouped_files_pairs[0]}_v_${grouped_files_pairs[2]}
    """
}
