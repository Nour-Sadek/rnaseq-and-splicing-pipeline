#!/usr/bin/env nextflow

/* Input csv file */
params.csv_reads = './test_data/test_csv.txt'

/* Output Directory */
// If user doesn't specify an output directory, let results be the default output directory
outputDir = params.outputDir ?: './results'

/* Parameters shared by at least two trimming tools */
params.base_quality_encoding = 'phred33'
params.min_len = 36

/* Trimmomatic parameters */
params.adapters_file = './adapters/TruSeq3-PE-2.fa'
params.seed_mismatches = 2
params.palindrome_clip_threshold = 30
params.simple_clip_threshold = 10
params.leading = 3
params.trailing = 3
params.window_size = 4
params.required_quality = 15

/* Bbduk parameters */
params.ktrim = 'r'
params.kmer_length = 23
params.min_kmer_size = 11
params.hdist = 1
params.tbo = true
params.tpe = true
tbo = params.tbo ? 'tbo' : ''
tpe = params.tpe ? 'tpe' : ''

/* Trim-Galore parameters */
params.quality = 20
params.stringency = 1

/* STAR parameters */
params.overhang = 151
params.genomeSAindexNbases = 12
params.filterMatch = 0.66

/* SALMON parameters */
params.kmer_size = 31

/* Default Parameters */
// If the user specifies a different parameter, it will be what the user specifies; these defaults only apply 
// if the user doesn't specify a value to the paramaters
params.genomeFastaFile = './genome_files/caenorhabditis_elegans.PRJNA13758.WBPS19.genomic.fa'
params.transcriptFastaFile = './genome_files/caenorhabditis_elegans.PRJNA13758.WBPS19.mRNA_transcripts.fa'
params.annotationsGTFFile = './genome_files/caenorhabditis_elegans.PRJNA13758.WBPS19.canonical_geneset.gtf'
params.type = 'paired-end'  // other option is 'single-end'
params.trimming = 'trimmomatic'  // other options are 'bbduk', 'trim_galore'
params.aligner = 'star'  // other options are 'hisat2', 'none'
params.quantifier = 'htseq-count'  // other options are 'featureCounts', 'htseq-count', 'salmon-alignment-mode', 'salmon-quasi-mapping-mode', 'kallisto', 'rsem'
params.splicingAnalyzer = 'majiq'  // other options are 'jum'

/* Included Modules */
include { FASTQC; FASTQC as FASTQCAFTERTRIMMING } from './modules/fastqc.nf'
include { TRIMMOMATIC; BBDUK; TRIM_GALORE } from './modules/trimming.nf'
include { STAR_REFERENCE_INDEX; HISAT2_REFERENCE_INDEX; HISAT2_IMPROVE_SPLICE_ALIGNMENT; SALMON_REFERENCE_INDEX } from './modules/reference_genome_index.nf'
include { STAR; HISAT2 } from './modules/aligning.nf'
include { SAM_TO_BAM; SORT_AND_INDEX_BAM } from './modules/samtools.nf'
include { HTSEQ_COUNT; FEATURE_COUNTS; SALMON_ALIGNMENT_MODE; SALMON_QUASI_MAPPING_MODE } from './modules/counting_reads.nf'

workflow {

    /* Read in the csv file */
    // Get the information in this format [sample_id, read_1, read_2]
    reads_channel = Channel.fromPath(params.csv_reads)
        .splitCsv(header:true)
        .map { row -> [row.sample_id, file(row.fastq_1), file(row.fastq_2)] }
    
    /* Run initial FastQC */
    fastqc_before_dir = outputDir + "/fastqc/before_trimming"
    FASTQC(reads_channel, fastqc_before_dir)

    /* Trim the reads */
    if (params.trimming == 'trimmomatic') {
        // Add all the trimmomatic arguments required into the <trimmomaticArgs> variable
        trimmomaticArgs = "${params.seed_mismatches}:${params.palindrome_clip_threshold}:${params.simple_clip_threshold} "
        trimmomaticArgs = trimmomaticArgs + "-${params.base_quality_encoding} LEADING:${params.leading} TRAILING:${params.trailing} "
        trimmomaticArgs = trimmomaticArgs + "SLIDINGWINDOW:${params.window_size}:${params.required_quality} MINLEN:${params.min_len}"
        
        // Run the TRIMMOMATIC process
        TRIMMOMATIC(reads_channel, file(params.adapters_file), outputDir, trimmomaticArgs)

        // Extract the sample_id, fwd_trimmed, and rev_trimmed outputs
        paired = TRIMMOMATIC.out.sample_id.combine(TRIMMOMATIC.out.fwd_trimmed)
        paired = paired.combine(TRIMMOMATIC.out.rev_trimmed)
        trimming_output_channel = paired.map { sample_id, fwd_trimmed, rev_trimmed -> [sample_id, fwd_trimmed, rev_trimmed] }
    } else if (params.trimming == 'bbduk') {
        // Add all the bbduk arguments required into the <bbdukArgs> variable
        bbdukArgs = "ktrim=${params.ktrim} k=${params.kmer_length} mink=${params.min_kmer_size} hdist=${params.hdist} "
        bbdukArgs = bbdukArgs + "minlen=${params.min_len} ${tbo} ${tpe}"

        // Run the BBDUK process
        BBDUK(reads_channel, outputDir, bbdukArgs)

        // Extract the sample_id, fwd_trimmed, and rev_trimmed outputs
        paired = BBDUK.out.sample_id.combine(BBDUK.out.fwd_trimmed)
        paired = paired.combine(BBDUK.out.rev_trimmed)
        trimming_output_channel = paired.map { sample_id, fwd_trimmed, rev_trimmed -> [sample_id, fwd_trimmed, rev_trimmed] }
    } else if (params.trimming == 'trim_galore') {
        // Add all the trim_galore arguments required into the <trimGaloreArgs> variable
        trimGaloreArgs = "--quality ${params.quality} --length ${params.min_len} --stringency ${params.stringency} --${params.base_quality_encoding}"

        // Run the TRIM_GALORE process
        TRIM_GALORE(reads_channel, outputDir, trimGaloreArgs)

        // Extract the sample_id, fwd_trimmed, and rev_trimmed outputs
        paired = TRIM_GALORE.out.sample_id.combine(TRIM_GALORE.out.fwd_trimmed)
        paired = paired.combine(TRIM_GALORE.out.rev_trimmed)
        trimming_output_channel = paired.map { sample_id, fwd_trimmed, rev_trimmed -> [sample_id, fwd_trimmed, rev_trimmed] }
    }

    /* Run fastqc after trimming */
    fastqc_after_dir = outputDir + "/fastqc/after_trimming"
    FASTQCAFTERTRIMMING(trimming_output_channel, fastqc_after_dir)

    /* Run the Alignment (if needed) */
    if (params.aligner == "star") {
        // Create the reference genome index
        STAR_REFERENCE_INDEX(outputDir, file(params.genomeFastaFile), file(params.annotationsGTFFile), params.overhang, params.genomeSAindexNbases)
        star_index_dir = outputDir + "/star/reference_index"

        // Run the STAR alignment process
        STAR(trimming_output_channel, outputDir, file(star_index_dir), params.filterMatch)

        // Extract the sample_id and bam_file outputs
        paired = STAR.out.sample_id.combine(STAR.out.bam_file)
        alignment_output_channel = paired.map { sample_id, bam_file -> [sample_id, bam_file] }
    } else if (params.aligner == 'hisat2') {
        // Create the reference genome index
        HISAT2_REFERENCE_INDEX(outputDir, file(params.genomeFastaFile))
        hisat2_index_dir = outputDir + "/hisat2/reference_index/genome"

        // Extract splice sites and exons for improved splicing alignment
        HISAT2_IMPROVE_SPLICE_ALIGNMENT(outputDir, file(params.annotationsGTFFile))

        // Run the HISAT2 alignment process
        HISAT2(trimming_output_channel, outputDir, file(hisat2_index_dir), HISAT2_IMPROVE_SPLICE_ALIGNMENT.out.splice_sites, HISAT2_IMPROVE_SPLICE_ALIGNMENT.out.exons)

        // Convert the sam_file to bam_file
        paired = HISAT2.out.sample_id.combine(HISAT2.out.sam_file)
        paired = paired.map { sample_id, sam_file -> [sample_id, sam_file] }
        SAM_TO_BAM(outputDir, paired)

        // Extract the sample_id and bam_file outputs
        paired = SAM_TO_BAM.out.sample_id.combine(SAM_TO_BAM.out.bam_file)
        alignment_output_channel = paired.map { sample_id, bam_file -> [sample_id, bam_file] }
    }

    /* Quantify the reads of the genes */
    if (params.aligner != 'none') {  // use quantifiers that need prior alignment
        /* Sort then index the bam files */
         SORT_AND_INDEX_BAM(alignment_output_channel, outputDir)
        paired = SORT_AND_INDEX_BAM.out.sample_id.combine(SORT_AND_INDEX_BAM.out.sorted_bam_file)
        sorted_bam_output_channel = paired.map { sample_id, sorted_bam_file -> [sample_id, sorted_bam_file] }

        /* Count number of reads for features (genes) using the sorted alignment bam files */
        if (params.quantifier == 'htseq-count') {
            // Run the HTSEQ_COUNT reads quantification process
            HTSEQ_COUNT(sorted_bam_output_channel, outputDir, file(params.annotationsGTFFile))
        } else if (params.quantifier == 'featureCounts') {
            // Run the FEATURE_COUNTS reads quantification process
            FEATURE_COUNTS(sorted_bam_output_channel, outputDir, file(params.annotationsGTFFile))
        } else if (params.quantifier == 'salmon-alignment-mode') {
            // Run the SALMON_ALIGNMENT_MODE reads quantification process
            SALMON_ALIGNMENT_MODE(sorted_bam_output_channel, outputDir, file(params.transcriptFastaFile))
        }
    } else {  // use quantifiers that don't need prior alignment
        if (params.quantifier == 'salmon-quasi-mapping-mode') {
            // Create the reference genome index
            SALMON_REFERENCE_INDEX(outputDir, file(params.transcriptFastaFile), params.kmer_size)
            salmon_index_dir = outputDir + "./salmon/reference_index"

            // Run the SALMON_QUASSI_MAPPING_MODE reads quantification process
            SALMON_QUASI_MAPPING_MODE(trimming_output_channel, outputDir, file(salmon_index_dir))
        } else if (params.quantifier == 'kallisto') {
            // Create the reference genome index

            // Run the KALLISTO reads quantification process
        }
    }
    






}