#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Specifying outputDir
def outputDir = params.outputDir ?: './results'

/* Included Modules */
include { FASTQC; FASTQC as FASTQCAFTERTRIMMING } from './modules/fastqc.nf'
include { BOWTIE2_CONTAMINANT_INDEX; BOWTIE2_REMOVE_CONTAMINANTS } from './modules/contaminant_rna_filtering.nf'
include { TRIMMOMATIC; BBDUK; TRIM_GALORE } from './modules/trimming.nf'
include { STAR_REFERENCE_INDEX; HISAT2_REFERENCE_INDEX; SALMON_REFERENCE_INDEX; KALLISTO_REFERENCE_INDEX; RSEM_REFERENCE_INDEX } from './modules/reference_genome_index.nf'
include { STAR; HISAT2 } from './modules/aligning.nf'
include { HTSEQ_COUNT; FEATURE_COUNTS; SALMON_QUASI_MAPPING_MODE; KALLISTO; RSEM } from './modules/counting_reads.nf'
include { SAM_TO_BAM; SORT_AND_INDEX_BAM } from './modules/samtools.nf'
include { MAJIQ_CONFIG; MAJIQ_BUILD; MAJIQ_PSI; MAJIQ_DELTA_PSI } from './modules/splicing/majiq.nf'
include { VOILA_PSI; VOILA_DELTA_PSI } from './modules/splicing/voila.nf'
include { rMATS_DIFFERENTIAL; rMATS_INDIVIDUAL } from './modules/splicing/rMATS.nf'
include { SUPPA2_GENERATE_EVENT_ANNOTATIONS; SUPPA2_CALCULATE_EVENTS_PSI; SUPPA2_SPLIT_FILES; SUPPA2_CALCULATE_EVENTS_DELTA_PSI } from './modules/splicing/suppa2.nf'
include { WHIPPET_INDEX; WHIPPET_QUANT; WHIPPET_DELTA } from './modules/splicing/whippet.nf'

workflow {

    /* Read in the csv file */
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
    if (params.run_fastqc_before_trimming) {
        fastqc_before_dir = outputDir + "/fastqc/before_trimming"
        FASTQC(reads_channel, fastqc_before_dir)
    }

    /* Clean contaminants */
    if (params.filter_contaminants && params.contaminants) {
        // Create the contaminants index
        contaminants = Channel.value(params.contaminants.collect { file(it) })  // e.g. [file('trna.fa'), file('ncrna.fa')]
        BOWTIE2_CONTAMINANT_INDEX(contaminants, outputDir)

        // Clean the fastq files from contaminants
        BOWTIE2_REMOVE_CONTAMINANTS(reads_channel, BOWTIE2_CONTAMINANT_INDEX.out.contaminants_index_files, outputDir)
        reads_channel = BOWTIE2_REMOVE_CONTAMINANTS.out.cleaned_samples
    }

    /* Trim the reads */
    if (params.trimming && params.trimming == 'trimmomatic') {

        // Run the TRIMMOMATIC process
        trimmomaticArgs = OrganizeArguments.makeTrimmomaticArgs(params.base_quality_encoding, params.adapters_file, params.seed_mismatches, params.palindrome_clip_threshold, params.simple_clip_threshold, 
                    params.min_adapter_length_palindrome, params.keepbothreads, params.window_size, params.required_quality, params.target_length, params.strictness, params.bases, 
                    params.min_count, params.max_count, params.leading, params.trailing, params.headcrop, params.tailcrop, params.crop, params.minlen, params.maxlen, params.avgqual)
        TRIMMOMATIC(reads_channel, outputDir, trimmomaticArgs)

        // Extract the sample_id, fwd_trimmed, and rev_trimmed outputs
        trimming_output_channel = TRIMMOMATIC.out.trimmed_samples
        
    } else if (params.trimming && params.trimming == 'bbduk') {

        // Run the BBDUK process
        bbdukArgs = OrganizeArguments.makeBbdukArgs(params.qin, params.reads, params.samplerate, params.k, params.rcomp, params.maskmiddle, params.minkmerhits, params.minkmerfraction, params.mincovfraction, 
                params.hammingdistance, params.qhdist, params.editdistance, params.hammingdistance2, params.qhdist2, params.editdistance2, params.forbidn, params.ktrim, params.ktrimtips, params.maskfullycovered, 
                params.mink, params.qtrim, params.trimq, params.minlength, params.minlengthfraction, params.minavgquality, params.minbasequality, params.maxns, params.minconsecutivebases, params.trimpad, 
                params.trimbyoverlap, params.strictoverlap, params.minoverlap, params.mininsert, params.trimpairsevenly, params.forcetrimleft, params.forcetrimright, params.forcetrimright2, params.forcetrimmod, 
                params.restrictleft, params.restrictright, params.mingc, params.maxgc, params.tossjunk)
        BBDUK(reads_channel, outputDir, bbdukArgs)

        // Extract the sample_id, fwd_trimmed, and rev_trimmed outputs
        trimming_output_channel = BBDUK.out.trimmed_samples

    } else if (params.trimming && params.trimming == 'trim_galore') {
        
        // Run the TRIM_GALORE process
        trimGaloreArgs = OrganizeArguments.makeTrimGaloreArgs(params.paired_end, params.quality, params.quality_encoding, params.adapter_sequence_1, params.adapter_sequence_2, params.specific_adapters, 
                params.max_length, params.stringency, params.error_rate, params.length, params.maxn, params.trim_n, params.trim_1, params.clip_R1, params.clip_R2, params.three_prime_clip_R1, params.three_prime_clip_R2, 
                params.nextseq_quality, params.hardtrim5, params.hardtrim3)
        TRIM_GALORE(reads_channel, outputDir, trimGaloreArgs)

        // Extract the sample_id, fwd_trimmed, and rev_trimmed outputs
        trimming_output_channel = TRIM_GALORE.out.trimmed_samples
    }

    /* Run fastqc after trimming */
    if (params.run_fastqc_after_trimming) {
        fastqc_after_dir = outputDir + "/fastqc/after_trimming"
        FASTQCAFTERTRIMMING(trimming_output_channel, fastqc_after_dir)
    }

    /* Run the Alignment (if needed) */
    if (params.aligner && params.aligner == "star") {

        // Create the reference genome index
        starReferenceIndexArgs = OrganizeArguments.makeStarReferenceIndexArgs(params.runRNGseed, params.genomeChrBinNbits, params.genomeSAindexNbases, params.genomeSAsparseD, params.genomeSuffixLengthMax, 
            params.sjdbGTFfeatureExon, params.sjdbGTFtagExonParentTranscript, params.sjdbGTFtagExonParentGene, params.sjdbGTFtagExonParentGeneName, params.sjdbGTFtagExonParentGeneType, params.sjdbOverhang, 
            params.sjdbScore, params.sjdbInsertSave)
        STAR_REFERENCE_INDEX(outputDir, file(params.genomeFastaFile), file(params.annotationsGTFFile), starReferenceIndexArgs)

        // Run the STAR alignment process
        starArgs = OrganizeArguments.makeStarArgs(params.paired_end, params.runRNGseed, params.readMapNumber, params.readMatesLengthsIn, params.readQualityScoreBase, params.clipAdapterType, params.clip3pNbases, 
            params.clip3pAdapterSeq, params.clip3pAdapterMMp, params.clip3pAfterAdapterNbases, params.clip5pNbases, params.outReadsUnmapped, params.outQSconversionAdd, params.outMultimapperOrder, params.outSAMmode, 
            params.outSAMstrandField, params.outSAMattributes, params.outSAMattrIHstart, params.outSAMunmapped, params.outSAMprimaryFlag, params.outSAMreadID, params.outSAMmapqUnique, params.outSAMflagOR, params.outSAMflagAND, 
            params.outSAMmultNmax, params.outSAMtlen, params.outBAMcompression, params.outFilterType, params.outFilterMultimapScoreRange, params.outFilterMultimapNmax, params.outFilterMismatchNmax, params.outFilterMismatchNoverLmax, 
            params.outFilterMismatchNoverReadLmax, params.outFilterScoreMin, params.outFilterScoreMinOverLread, params.outFilterMatchNmin, params.outFilterMatchNminOverLread, params.outFilterIntronMotifs, params.outFilterIntronStrands, 
            params.outSJfilterReads, params.outSJfilterOverhangMin, params.outSJfilterCountUniqueMin, params.outSJfilterCountTotalMin, params.outSJfilterDistToOtherSJmin, params.outSJfilterIntronMaxVsReadN, params.scoreGap, 
            params.scoreGapNoncan, params.scoreGapGCAG, params.scoreGapATAC, params.scoreGenomicLengthLog2scale, params.scoreDelOpen, params.scoreDelBase, params.scoreInsOpen, params.scoreInsBase, params.scoreStitchSJshift, 
            params.seedSearchStartLmax, params.seedSearchStartLmaxOverLread, params.seedSearchLmax, params.seedMultimapNmax, params.seedPerReadNmax, params.seedPerWindowNmax, params.seedNoneLociPerWindow, params.seedSplitMin, 
            params.seedMapMin, params.alignIntronMin, params.alignIntronMax, params.alignMatesGapMax, params.alignSJoverhangMin, params.alignSJstitchMismatchNmax, params.alignSJDBoverhangMin, params.alignSplicedMateMapLmin, 
            params.alignSplicedMateMapLminOverLmate, params.alignWindowsPerReadNmax, params.alignTranscriptsPerWindowNmax, params.alignTranscriptsPerReadNmax, params.alignEndsType, params.alignEndsProtrude, params.alignSoftClipAtReferenceEnds, 
            params.alignInsertionFlush, params.peOverlapNbasesMin, params.peOverlapMMp, params.winAnchorMultimapNmax, params.winBinNbits, params.winAnchorDistNbins, params.winFlankNbins, params.winReadCoverageRelativeMin, 
            params.winReadCoverageBasesMin, params.quantMode, params.quantTranscriptomeBAMcompression, params.quantTranscriptomeSAMoutput)
        STAR(trimming_output_channel, outputDir, STAR_REFERENCE_INDEX.out.reference_index, starArgs)

        // Extract the sample_id and bam_file outputs
        alignment_output_channel = STAR.out.alignment_output

    } else if (params.aligner && params.aligner == 'hisat2') {

        // Define the basename of the reference index
        def hisat2_index_prefix = params.hisat2_index_prefix ?: 'genome'

        // Create the reference genome index
        hisat2ReferenceIndexArgs = OrganizeArguments.makeHisat2ReferenceIndexArgs(params.large_index, params.noauto, params.bmax, params.bmaxdivn, params.dcv, params.nodc, params.noref, params.justref, params.offrate, 
            params.ftabchars, params.localoffrate, params.localftabchars, params.seed, params.cutoff)
        HISAT2_REFERENCE_INDEX(outputDir, hisat2_index_prefix, file(params.genomeFastaFile), hisat2ReferenceIndexArgs)

        // Run the HISAT2 alignment process
        hisat2Args = OrganizeArguments.makeHisat2Args(params.paired_end, params.skip, params.upto, params.trim5, params.trim3, params.phred_quality, params.solexa_quals, params.int_quals, params.n_ceil_func, 
            params.ignore_quals, params.norc, params.nofw, params.mismatch_penalties, params.soft_clipping, params.no_softclip, params.n_penalty, params.read_gap_penalty, params.reference_gap_penalty, 
            params.score_min_func, params.pen_cansplice, params.pen_noncansplice, params.pen_canintronlen, params.pen_noncanintronlen, params.min_intronlen, params.max_intronlen, params.no_temp_splicesite, 
            params.no_spliced_alignment, params.rna_strandness, params.transcriptome_mapping_only, params.downstream_transcriptome_assembly, params.dta_cufflinks, params.avoid_pseudogene, params.no_templatelen_adjustment, 
            params.num_alignments_per_read, params.max_seeds, params.report_all_alignments, params.report_secondary_alignments, params.min_fragment_length, params.max_fragment_length, params.mate_orientations, 
            params.no_mixed, params.no_discordant, params.index_offrate, params.reorder, params.rng_seed, params.non_deterministic)
        HISAT2(trimming_output_channel, outputDir, hisat2_index_prefix, HISAT2_REFERENCE_INDEX.out.hisat2_index_files, hisat2Args)

        // Convert the sam_file to bam_file
        SAM_TO_BAM(outputDir, HISAT2.out.alignment_output)

        // Extract the sample_id and bam_file outputs
        alignment_output_channel = SAM_TO_BAM.out.bam_output
    }

    /* Quantify the reads of the genes */
    if (params.aligner && params.aligner != 'none') {  // use quantifiers that need prior alignment
        /* Sort then index the bam files */
        SORT_AND_INDEX_BAM(alignment_output_channel, outputDir)
        sorted_bam_output_channel = SORT_AND_INDEX_BAM.out.sorted_bam_output

        /* Count number of reads for features (genes) using the sorted alignment bam files */
        if (params.quantifier && params.quantifier == 'htseq-count') {
            
            // Run the HTSEQ_COUNT reads quantification process
            htseqCountArgs = OrganizeArguments.makeHtseqCountArgs(params.max_reads_in_buffer, params.stranded, params.minaqual, params.feature_type, params.id_attribute, params.additional_attributes, params.mode, params.nonunique_mode, 
                params.secondary_alignments, params.supplementary_alignments, params.add_chromosome_info)
            HTSEQ_COUNT(sorted_bam_output_channel, outputDir, file(params.annotationsGTFFile), htseqCountArgs)

        } else if (params.quantifier && params.quantifier == 'featureCounts') {
            
            // Run the FEATURE_COUNTS reads quantification process
            featureCountsArgs = OrganizeArguments.makeFeatureCountsArgs(params.paired_end, params.requireBothEndsMapped, params.countChimericFragments, params.checkFragLength, params.countReadPairs, params.autosort, params.minFragLength, 
                params.maxfragLength, params.useMetaFeatures, params.isGTFAnnotationFile, params.attrType_GTF, params.juncCounts, params.isLongRead, params.countMultiMappingReads, params.allowMultiOverlap, params.minMQS, params.isStrandSpecific, 
                params.featureType_GTF, params.byReadGroup, params.attrType_GTF_extra, params.fraction, params.fracOverlap, params.fracOverlapFeature, params.ignoreDup, params.largestOverlap, params.minOverlap, params.nonOverlap, 
                params.nonOverlapFeature, params.nonSplitOnly, params.primaryOnly, params.read2pos, params.readExtension3, params.readExtension5, params.readShiftSize, params.readShiftType, params.splitOnly)
            FEATURE_COUNTS(sorted_bam_output_channel, outputDir, file(params.annotationsGTFFile), featureCountsArgs)
        }

        /* Perform splicing alignment analysis */
        if (params.splicingAnalyzer && params.splicingAnalyzer == 'majiq') {
            
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

        } else if (params.splicingAnalyzer && params.splicingAnalyzer == 'rMats') {

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
        if (params.quantifier && params.quantifier == 'salmon-quasi-mapping-mode') {
            // Create the reference genome index
            SALMON_REFERENCE_INDEX(outputDir, file(params.transcriptFastaFile), params.kmer_size)

            // Run the SALMON_QUASSI_MAPPING_MODE reads quantification process
            SALMON_QUASI_MAPPING_MODE(trimming_output_channel, outputDir, SALMON_REFERENCE_INDEX.out.reference_index)
            tpm_column = 4
            
            quantifier_output_channel = SALMON_QUASI_MAPPING_MODE.out.quants_file
        } else if (params.quantifier && params.quantifier == 'kallisto') {
            // Create the reference genome index
            KALLISTO_REFERENCE_INDEX(outputDir, file(params.transcriptFastaFile))

            // Run the KALLISTO reads quantification process
            KALLISTO(trimming_output_channel, outputDir, KALLISTO_REFERENCE_INDEX.out.kallisto_reference_index, params.num_bootstrap_samples)
            tpm_column = 5

            quantifier_output_channel = KALLISTO.out.quants_file
        } else if (params.quantifier && params.quantifier == 'rsem') {
            // Create the reference genome index
            RSEM_REFERENCE_INDEX(outputDir, params.rsem_index_prefix, file(params.genomeFastaFile), file(params.annotationsGTFFile), params.overhang, params.genomeSAindexNbases)

            // Run the RSEM alignment and quantification process
            RSEM(trimming_output_channel, outputDir, params.rsem_index_prefix, RSEM_REFERENCE_INDEX.out.rsem_index_files)
        }

        if (params.splicingAnalyzer && params.splicingAnalyzer == 'suppa2') {
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
                r_file = file("./modules/splicing/split_file.R")
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
        } else if (params.splicingAnalyzer && params.splicingAnalyzer == 'whippet') {
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
