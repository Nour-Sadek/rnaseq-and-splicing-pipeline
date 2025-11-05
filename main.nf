#!/usr/bin/env nextflow

/* Input csv file */
params.csv_reads = './test_data_big/test_csv.txt'

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
params.overhang = 101
params.genomeSAindexNbases = 12
params.filterMatch = 0.3

/* HISAT2 parameters */
params.hisat2_index_prefix = 'genome'

/* SALMON parameters */
params.kmer_size = 31

/* KALLISTO parameters */
params.num_bootstrap_samples = 0

/* RSEM parameters */
params.rsem_index_prefix = 'rsem_reference'

/* General Splicing parameters */
params.individualSplicingAnalysis = true
params.differentialSplicingAnalysis = true

/* rMats parameters */
params.read_length = 100


/* Default Parameters */
// If the user specifies a different parameter, it will be what the user specifies; these defaults only apply 
// if the user doesn't specify a value to the paramaters
params.genomeFastaFile = './genome_files/GCF_000002985.6_WBcel235_genomic.fna'
params.transcriptFastaFile = './genome_files/rna.fna'
params.annotationsGTFFile = './genome_files/genomic.gtf'
params.annotationsGFF3File = './genome_files/genomic.gff'
params.trimming = 'trimmomatic'  // other options are 'bbduk', 'trim_galore'
params.aligner = 'star'  // other options are 'hisat2', 'minimap2', 'none'
params.quantifier = 'htseq-count'  // other options are 'featureCounts', 'htseq-count', 'salmon-quasi-mapping-mode', 'kallisto', 'rsem'
params.splicingAnalyzer = 'majiq'  // other options are 'rMats', 'suppa2', 'whippet', 'majiq'

/* Included Modules */
include { FASTQC; FASTQC as FASTQCAFTERTRIMMING } from './modules/fastqc.nf'
include { TRIMMOMATIC; BBDUK; TRIM_GALORE } from './modules/trimming.nf'
include { STAR_REFERENCE_INDEX; HISAT2_REFERENCE_INDEX; MINIMAP2_REFERENCE_INDEX; SALMON_REFERENCE_INDEX; KALLISTO_REFERENCE_INDEX; RSEM_REFERENCE_INDEX } from './modules/reference_genome_index.nf'
include { STAR; HISAT2; MINIMAP2 } from './modules/aligning.nf'
include { SAM_TO_BAM; SORT_AND_INDEX_BAM } from './modules/samtools.nf'
include { GFFREAD } from './modules/gff_utilities.nf'
include { HTSEQ_COUNT; FEATURE_COUNTS; SALMON_ALIGNMENT_MODE; SALMON_QUASI_MAPPING_MODE; KALLISTO; RSEM } from './modules/counting_reads.nf'
include { MAJIQ_CONFIG; MAJIQ_BUILD; MAJIQ_PSI; MAJIQ_DELTA_PSI; VOILA_PSI; rMATS_DIFFERENTIAL; rMATS_INDIVIDUAL } from './modules/splicing.nf'

workflow {

    /* Read in the csv file */
    // Get the information in this format [sample_id, sample_group, read_1, read_2]
    reads_channel = Channel.fromPath(params.csv_reads)
        .splitCsv(header:true)
        .map { row -> [row.sample_id, row.sample_group, file(row.fastq_1), file(row.fastq_2)] }
    
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
        trimming_output_channel = TRIMMOMATIC.out.trimmed_samples
    } else if (params.trimming == 'bbduk') {
        // Add all the bbduk arguments required into the <bbdukArgs> variable
        bbdukArgs = "ktrim=${params.ktrim} k=${params.kmer_length} mink=${params.min_kmer_size} hdist=${params.hdist} "
        bbdukArgs = bbdukArgs + "minlen=${params.min_len} ${tbo} ${tpe}"

        // Run the BBDUK process
        BBDUK(reads_channel, outputDir, bbdukArgs)

        // Extract the sample_id, fwd_trimmed, and rev_trimmed outputs
        trimming_output_channel = BBDUK.out.trimmed_samples
    } else if (params.trimming == 'trim_galore') {
        // Add all the trim_galore arguments required into the <trimGaloreArgs> variable
        trimGaloreArgs = "--quality ${params.quality} --length ${params.min_len} --stringency ${params.stringency} --${params.base_quality_encoding}"

        // Run the TRIM_GALORE process
        TRIM_GALORE(reads_channel, outputDir, trimGaloreArgs)

        // Extract the sample_id, fwd_trimmed, and rev_trimmed outputs
        trimming_output_channel = TRIM_GALORE.out.trimmed_samples
    }

    /* Run fastqc after trimming */
    fastqc_after_dir = outputDir + "/fastqc/after_trimming"
    FASTQCAFTERTRIMMING(trimming_output_channel, fastqc_after_dir)

    /* Run the Alignment (if needed) */
    if (params.aligner == "star") {
        // Create the reference genome index
        STAR_REFERENCE_INDEX(outputDir, file(params.genomeFastaFile), file(params.annotationsGTFFile), params.overhang, params.genomeSAindexNbases)

        // Run the STAR alignment process
        STAR(trimming_output_channel, outputDir, STAR_REFERENCE_INDEX.out.reference_index, params.filterMatch)

        // Extract the sample_id and bam_file outputs
        alignment_output_channel = STAR.out.alignment_output
    } else if (params.aligner == 'hisat2') {
        // Create the reference genome index
        HISAT2_REFERENCE_INDEX(outputDir, params.hisat2_index_prefix, file(params.genomeFastaFile))

        // Run the HISAT2 alignment process
        HISAT2(trimming_output_channel, outputDir, params.hisat2_index_prefix, HISAT2_REFERENCE_INDEX.out.hisat2_index_files)

        // Convert the sam_file to bam_file
        SAM_TO_BAM(outputDir, HISAT2.out.alignment_output)

        // Extract the sample_id and bam_file outputs
        alignment_output_channel = SAM_TO_BAM.out.bam_output
    }

    /* Quantify the reads of the genes */
    if (params.aligner != 'none') {  // use quantifiers that need prior alignment
        /* Sort then index the bam files */
        SORT_AND_INDEX_BAM(alignment_output_channel, outputDir)
        sorted_bam_output_channel = SORT_AND_INDEX_BAM.out.sorted_bam_output

        /* Count number of reads for features (genes) using the sorted alignment bam files */
        if (params.quantifier == 'htseq-count') {
            // Run the HTSEQ_COUNT reads quantification process
            HTSEQ_COUNT(sorted_bam_output_channel, outputDir, file(params.annotationsGTFFile))
        } else if (params.quantifier == 'featureCounts') {
            // Run the FEATURE_COUNTS reads quantification process
            FEATURE_COUNTS(sorted_bam_output_channel, outputDir, file(params.annotationsGTFFile))
        }

        /* Perform splicing alignment analysis */
        if (params.splicingAnalyzer == 'majiq') {
            
            // Build the config file
            all_sample_bams = SORT_AND_INDEX_BAM.out.sorted_bam_output
                    .groupTuple(by: 1)  
                    .map { samples, group, reads -> [group, reads.simpleName] }.collect(flat: false, sort: true)
            
            bam_dirs = SORT_AND_INDEX_BAM.out.sorted_bam_output
                .groupTuple(by: 1)  
                .map { samples, group, reads -> [group, reads] }.collect(flat: false, sort: true)
                .collectMany { it[1] }  // get only the bam paths     
                .collect { it.parent }  // get the parent work folder for each bam file
                .unique()
            
            MAJIQ_CONFIG(all_sample_bams, bam_dirs, outputDir)

            // Run majiq build
            MAJIQ_BUILD(MAJIQ_CONFIG.out.config, file(params.annotationsGFF3File), outputDir)

            if (params.individualSplicingAnalysis) {
                // Create a channel that would return values such as [sample_group, [files names of replicates]]
                groups_file_names = SORT_AND_INDEX_BAM.out.sorted_bam_output
                    .groupTuple(by: 1)  
                    .map { samples, group, reads -> tuple(group, reads.simpleName) }
                
                // Run majiq psi
                MAJIQ_PSI(groups_file_names, MAJIQ_BUILD.out.samples_splice_graphs, outputDir)

                // Get the work folder for the psi files
                majiq_psi_parent_folder = MAJIQ_PSI.out.majiq_tsv_file.map { it.parent }.view()

                // Run voila psi
                VOILA_PSI(MAJIQ_PSI.out.sample_group, majiq_psi_parent_folder, MAJIQ_BUILD.out.splicegraph_file, outputDir)
            }

            if (params.differentialSplicingAnalysis) {
                // Return a channel were each output is of the format:
                // [sample_group_1, [sample_group_1's replicates file names], sample_group_2, [sample_group_2's replicates file names]]
                grouped_files_pairs = SORT_AND_INDEX_BAM.out.sorted_bam_output
                    .groupTuple(by: 1, sort: true)  
                    .map { samples, group, reads -> [group, reads.simpleName] }
                    .toList()
                    .flatMap { grouped_list ->
                        def pairs = []
                        for (i in 0..<grouped_list.size()) {
                            for (j in i+1..<grouped_list.size()) {
                                pairs << [grouped_list[i][0], grouped_list[i][1], grouped_list[j][0], grouped_list[j][1]]
                            }
                        }
                        return pairs
                    }
                
                // Run differential splicing analysis on each possible pair of samples
                MAJIQ_DELTA_PSI(grouped_files_pairs, MAJIQ_BUILD.out.samples_splice_graphs, outputDir)
            }

        } else if (params.splicingAnalyzer == 'rMats') {

            if (params.individualSplicingAnalysis) {
                // Return a channel where each output is of the format:
                // [sample_group, [bam_files of replicates]]
                sample_bams = SORT_AND_INDEX_BAM.out.sorted_bam_output
                    .groupTuple(by: 1)  
                    .map { samples, group, reads -> [group, reads] }
                
                // Run splicing analysis on each sample group
                rMATS_INDIVIDUAL(sample_bams, params.read_length, file(params.annotationsGTFFile), outputDir)
            }

            if (params.differentialSplicingAnalysis) {
                // Return a channel were each output is of the format:
                // [[sample_group_1, [sample_group_1's replicates bam_files]], [sample_group_2, [sample_group_2's replicates bam_files]]]
                grouped_bams_pairs = SORT_AND_INDEX_BAM.out.sorted_bam_output
                    .groupTuple(by: 1, sort: true)  
                    .map { samples, group, reads -> [group, reads] } // each channel would be in this form [sample_group, [bam_files of replicates]]
                    .toList()
                    .flatMap { grouped_list ->
                        def pairs = []
                        for (i in 0..<grouped_list.size()) {
                            for (j in i+1..<grouped_list.size()) {
                                pairs << [grouped_list[i][0], grouped_list[i][1], grouped_list[j][0], grouped_list[j][1]]
                            }
                        }
                        return pairs
                    }
            
                // Run differential splicing analysis on each possible pair of samples
                rMATS_DIFFERENTIAL(grouped_bams_pairs, params.read_length, file(params.annotationsGTFFile), outputDir)
            }


        } else if (params.splicingAnalyzer == 'suppa2') {

        } else if (params.splicingAnalyzer == 'whippet') {

        }

    } else {  // use quantifiers that don't need prior alignment
        if (params.quantifier == 'salmon-quasi-mapping-mode') {
            // Create the reference genome index
            SALMON_REFERENCE_INDEX(outputDir, file(params.transcriptFastaFile), params.kmer_size)

            // Run the SALMON_QUASSI_MAPPING_MODE reads quantification process
            SALMON_QUASI_MAPPING_MODE(trimming_output_channel, outputDir, SALMON_REFERENCE_INDEX.out.reference_index)
        } else if (params.quantifier == 'kallisto') {
            // Create the reference genome index
            KALLISTO_REFERENCE_INDEX(outputDir, file(params.transcriptFastaFile))

            // Run the KALLISTO reads quantification process
            KALLISTO(trimming_output_channel, outputDir, KALLISTO_REFERENCE_INDEX.out.kallisto_reference_index, params.num_bootstrap_samples)
        } else if (params.quantifier == 'rsem') {
            // Create the reference genome index
            RSEM_REFERENCE_INDEX(outputDir, params.rsem_index_prefix, file(params.genomeFastaFile), file(params.annotationsGTFFile))

            // Run the RSEM alignment and quantification process
            RSEM(trimming_output_channel, outputDir, params.rsem_index_prefix, RSEM_REFERENCE_INDEX.out.rsem_index_files)
        }
    }
    






}