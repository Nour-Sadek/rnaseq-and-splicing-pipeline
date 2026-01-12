#!/usr/bin/env nextflow

// Specifying outputDir and specific bbduk parameters based on the parameters from the yaml file
outputDir = params.outputDir ?: './results'

// Specify whether reasd are paired or single end
params.paired_end = false

/* Included Modules */
include { FASTQC; FASTQC as FASTQCAFTERTRIMMING } from './modules/fastqc.nf'
include { TRIMMOMATIC; BBDUK; TRIM_GALORE } from './modules/trimming.nf'
include { STAR_REFERENCE_INDEX; HISAT2_REFERENCE_INDEX; SALMON_REFERENCE_INDEX; KALLISTO_REFERENCE_INDEX; RSEM_REFERENCE_INDEX } from './modules/reference_genome_index.nf'
include { STAR; HISAT2 } from './modules/aligning.nf'
include { HTSEQ_COUNT; FEATURE_COUNTS; SALMON_QUASI_MAPPING_MODE; KALLISTO; RSEM } from './modules/counting_reads.nf'
include { SAM_TO_BAM; SORT_AND_INDEX_BAM } from './modules/samtools.nf'
// processes for splicing
include { MAJIQ_CONFIG; MAJIQ_BUILD; MAJIQ_PSI; MAJIQ_DELTA_PSI } from './modules/splicing/majiq.nf'
include { VOILA_PSI; VOILA_DELTA_PSI } from './modules/splicing/voila.nf'
include { rMATS_DIFFERENTIAL; rMATS_INDIVIDUAL } from './modules/splicing/rMATS.nf'
include { SUPPA2_GENERATE_EVENT_ANNOTATIONS; SUPPA2_CALCULATE_EVENTS_PSI; SUPPA2_SPLIT_FILES; SUPPA2_CALCULATE_EVENTS_DELTA_PSI } from './modules/splicing/suppa2.nf'
include { WHIPPET_INDEX; WHIPPET_QUANT; WHIPPET_DELTA } from './modules/splicing/whippet.nf'

workflow {

    /* Read in the csv file */
    // Get the information in this format [sample_id, sample_group, read_1, read_2] or [sample_id, sample_group, read]
    if (params.paired_end) {
        reads_channel = Channel.fromPath(params.csv_reads)
        .splitCsv(header:true)
        .map { row -> [row.sample_id, row.sample_group, [file(row.fastq_1), file(row.fastq_2)]] }
    } else {
        reads_channel = Channel.fromPath(params.csv_reads)
        .splitCsv(header:true)
        .map { row -> [row.sample_id, row.sample_group, [file(row.fastq_1)]] }
    }
    
    /* Run initial FastQC */
    // fastqc_before_dir = outputDir + "/fastqc/before_trimming"
    // FASTQC(reads_channel, fastqc_before_dir)

    /* Trim the reads */
    if (params.trimming == 'trimmomatic') {
        // Add all the trimmomatic arguments required into the <trimmomaticArgs> variable
        trimmomaticArgs = ["-${params.base_quality_encoding}"]

        // Specifying the parameters for ILLUMINACLIP
        if (params.seed_mismatches != "none" && params.palindrome_clip_threshold != "none" && params.simple_clip_threshold != "none") {
            if (params.min_adapter_length_palindrome != "none" && params.keepbothreads != "none") {
                trimmomaticArgs << "ILLUMINACLIP:${params.adapters_file}:${params.seed_mismatches}:${params.palindrome_clip_threshold}:${params.simple_clip_threshold}:${params.min_adapter_length_palindrome}:${params.keepbothreads}"
            } else {
                trimmomaticArgs << "ILLUMINACLIP:${params.adapters_file}:${params.seed_mismatches}:${params.palindrome_clip_threshold}:${params.simple_clip_threshold}"
            }
            // /opt/conda/share/trimmomatic-0.40-0/adapters/
        }

        // Specifying the parameters for SLIDINGWINDOW
        if (params.windowsize != "none" && params.required_quality != "none") {
            trimmomaticArgs << "SLIDINGWINDOW:${params.window_size}:${params.required_quality}"
        }

        // Specifying the parameters for MAXINFO
        if (params.target_length != "none" && params.strictness != "none") {
            trimmomaticArgs << "MAXINFO:${params.target_length}:${params.strictness}"
        }

        // Specifying the parameters for BASECOUNT
        if (params.bases != "none" && params.min_count != "none" && params.max_count != "none") {
            trimmomaticArgs << "BASECOUNT:${params.bases}:${params.min_count}:${params.max_count}"
        }

        // Specifying the single parameters
        if (params.leading != "none") trimmomaticArgs << "LEADING:${params.leading}"
        if (params.trailing != "none") trimmomaticArgs << "TRAILING:${params.trailing}"
        if (params.headcrop != "none") trimmomaticArgs << "HEADCROP:${params.headcrop}"
        if (params.tailcrop != "none") trimmomaticArgs << "TAILCROP:${params.tailcrop}"
        if (params.crop != "none") trimmomaticArgs << "CROP:${params.crop}"
        if (params.minlen != "none") trimmomaticArgs << "MINLEN:${params.minlen}"
        if (params.maxlen != "none") trimmomaticArgs << "MAXLEN:${params.maxlen}"
        if (params.avgqual != "none") trimmomaticArgs << "AVGQUAL:${params.avgqual}"
        trimmomaticArgs = trimmomaticArgs.join(" ")

        // Run the TRIMMOMATIC process
        TRIMMOMATIC(reads_channel, outputDir, trimmomaticArgs)

        // Extract the sample_id, fwd_trimmed, and rev_trimmed outputs
        trimming_output_channel = TRIMMOMATIC.out.trimmed_samples
        
    } else if (params.trimming == 'bbduk') {
        // Add all the bbduk arguments required into the <bbdukArgs> variable
        bbdukArgs = [
            "qin=${params.qin}", "reads=${params.reads}", "samplerate=${params.samplerate}", "k=${params.k}", "rcomp=${params.rcomp}",
            "maskmiddle=${params.maskmiddle}", "minkmerhits=${params.minkmerhits}", "minkmerfraction=${params.minkmerfraction}",
            "mincovfraction=${params.mincovfraction}", "hammingdistance=${params.hammingdistance}", "editdistance=${params.editdistance}", 
            "hammingdistance2=${params.hammingdistance2}", "editdistance2=${params.editdistance2}", "forbidn=${params.forbidn}", 
            "ktrim=${params.ktrim}", "maskfullycovered=${params.maskfullycovered}", "mink=${params.mink}", "qtrim=${params.qtrim}", 
            "trimq=${params.trimq}", "minlength=${params.minlength}", "minlengthfraction=${params.minlengthfraction}", "minavgquality=${params.minavgquality}", 
            "minbasequality=${params.minbasequality}", "maxns=${params.maxns}", "minconsecutivebases=${params.minconsecutivebases}", 
            "trimpad=${params.trimpad}", "trimbyoverlap=${params.trimbyoverlap}", "strictoverlap=${params.strictoverlap}", "minoverlap=${params.minoverlap}", 
            "mininsert=${params.mininsert}", "trimpairsevenly=${params.trimpairsevenly}", "forcetrimleft=${params.forcetrimleft}",  
            "forcetrimright=${params.forcetrimright}", "forcetrimright2=${params.forcetrimright2}", "forcetrimmod=${params.forcetrimmod}",  
            "restrictleft=${params.restrictleft}", "restrictright=${params.restrictright}", "mingc=${params.mingc}", "maxgc=${params.maxgc}", 
            "tossjunk=${params.tossjunk}"
        ]

        if (params.mink == 0) bbdukArgs.remove("mink=${params.mink}")
        bbdukArgs = bbdukArgs.join(" ")

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
    //fastqc_after_dir = outputDir + "/fastqc/after_trimming"
    //FASTQCAFTERTRIMMING(trimming_output_channel, fastqc_after_dir)

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
                majiq_psi_parent_folder = MAJIQ_PSI.out.majiq_tsv_file.map { it.parent }

                // Run voila modulize psi
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

                // Get the work folder for the delta psi files
                majiq_delta_psi_parent_folder = MAJIQ_DELTA_PSI.out.majiq_tsv_file.map { it.parent }

                // Run voila modulize delta psi
                VOILA_DELTA_PSI(MAJIQ_DELTA_PSI.out.paired_samples_name, majiq_delta_psi_parent_folder, MAJIQ_BUILD.out.splicegraph_file, outputDir)
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

        }

    } else {  // use quantifiers that don't need prior alignment
        if (params.quantifier == 'salmon-quasi-mapping-mode') {
            // Create the reference genome index
            SALMON_REFERENCE_INDEX(outputDir, file(params.transcriptFastaFile), params.kmer_size)

            // Run the SALMON_QUASSI_MAPPING_MODE reads quantification process
            SALMON_QUASI_MAPPING_MODE(trimming_output_channel, outputDir, SALMON_REFERENCE_INDEX.out.reference_index)
            tpm_column = 4
            
            quantifier_output_channel = SALMON_QUASI_MAPPING_MODE.out.quants_file
        } else if (params.quantifier == 'kallisto') {
            // Create the reference genome index
            KALLISTO_REFERENCE_INDEX(outputDir, file(params.transcriptFastaFile))

            // Run the KALLISTO reads quantification process
            KALLISTO(trimming_output_channel, outputDir, KALLISTO_REFERENCE_INDEX.out.kallisto_reference_index, params.num_bootstrap_samples)
            tpm_column = 5

            quantifier_output_channel = KALLISTO.out.quants_file
        } else if (params.quantifier == 'rsem') {
            // Create the reference genome index
            RSEM_REFERENCE_INDEX(outputDir, params.rsem_index_prefix, file(params.genomeFastaFile), file(params.annotationsGTFFile), params.overhang, params.genomeSAindexNbases)

            // Run the RSEM alignment and quantification process
            RSEM(trimming_output_channel, outputDir, params.rsem_index_prefix, RSEM_REFERENCE_INDEX.out.rsem_index_files)
        }

        if (params.splicingAnalyzer == 'suppa2') {
            // Generate the event annotations (ioe) file
            SUPPA2_GENERATE_EVENT_ANNOTATIONS(file(params.annotationsGTFFile), outputDir)

            if (params.individualSplicingAnalysis || params.differentialSplicingAnalysis) {
                all_salmon_samples = quantifier_output_channel.collect(flat: false)
                all_sample_ids = all_salmon_samples.map { it*.get(0) }
                all_sample_quants = all_salmon_samples.map { it*.get(2) }

                SUPPA2_CALCULATE_EVENTS_PSI(all_sample_ids, all_sample_quants, SUPPA2_GENERATE_EVENT_ANNOTATIONS.out.ioe_file, tpm_column, outputDir)
            }

            if (params.differentialSplicingAnalysis) {
                // Create a channel that would return values such as [sample_group, [files names of replicates]]
                groups_file_names_quant = quantifier_output_channel
                    .groupTuple(by: 1)  
                    .map { sample_id, group, quant -> tuple(group, sample_id) }
                
                // Seperate the psi and TPM files for the conditions
                r_file = file("./modules/split_file.R")
                SUPPA2_SPLIT_FILES(r_file, groups_file_names_quant, SUPPA2_CALCULATE_EVENTS_PSI.out.tpm_file, SUPPA2_CALCULATE_EVENTS_PSI.out.psi_file, outputDir)

                // Return a channel were each output is of the format:
                // [sample_group_1, sample_group_1_psi_file, sample_group_1_tpm_file, sample_group_2, sample_group_2_psi_file, sample_group_2_tpm_file]
                all_suppa2_files = SUPPA2_SPLIT_FILES.out.individual_sample_files.collect(flat: false)
                paired_suppa2 = all_suppa2_files
                    .flatMap { grouped_list ->
                        def pairs = []
                        for (i in 0..<grouped_list.size()) {
                            for (j in i+1..<grouped_list.size()) {
                                pairs << [grouped_list[i][0], grouped_list[i][1], grouped_list[i][2], grouped_list[j][0], grouped_list[j][1], grouped_list[j][2]]
                            }
                        }
                        return pairs
                    }
                
                // Perform differential splicing analysis
                SUPPA2_CALCULATE_EVENTS_DELTA_PSI(paired_suppa2, SUPPA2_GENERATE_EVENT_ANNOTATIONS.out.ioe_file, outputDir)
                
            }
        } else if (params.splicingAnalyzer == 'whippet') {
            // Build the index for Whippet
            WHIPPET_INDEX(file(params.genomeFastaFile), file(params.annotationsGTFFile), outputDir)

            if (params.individualSplicingAnalysis || params.differentialSplicingAnalysis) {
                WHIPPET_QUANT(trimming_output_channel, WHIPPET_INDEX.out.whippet_index, outputDir)
            }

            if (params.differentialSplicingAnalysis) {
                // Get the pairs of psi files between each sample group
                // The output of the channel would be [group_1, [group_1's replicates psi files], group_2 [group_2's replicates psu files]]
                grouped_files_pairs = WHIPPET_QUANT.out.sample_psi_file
                    .groupTuple(by: 1, sort: true)
                    .map { sample_id, group, psi_file -> [group, psi_file] }
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
                
                // Run the delta psi process
                all_samples_psi_files = WHIPPET_QUANT.out.sample_psi_file 
                    .map { sample_id, group, psi_file -> psi_file }.collect()
                WHIPPET_DELTA(grouped_files_pairs, all_samples_psi_files, outputDir)
            }
        }
    }

}
