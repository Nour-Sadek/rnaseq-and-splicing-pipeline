class OrganizeArguments {

    /* Making arguments for the trimming tools */
    static String makeTrimmomaticArgs(base_quality_encoding, adapters_file, seed_mismatches, palindrome_clip_threshold, simple_clip_threshold, 
    min_adapter_length_palindrome, keepbothreads, window_size, required_quality, target_length, strictness, bases, 
    min_count, max_count, leading, trailing, headcrop, tailcrop, crop, minlen, maxlen, avgqual) {

        // Add all the trimmomatic arguments required into the <trimmomaticArgs> variable
        def trimmomaticArgs = []
        if (base_quality_encoding) trimmomaticArgs << "-${base_quality_encoding}"

        // Specifying the parameters for ILLUMINACLIP
        if (seed_mismatches && palindrome_clip_threshold && simple_clip_threshold) {
            if (min_adapter_length_palindrome && keepbothreads) {
                trimmomaticArgs << "ILLUMINACLIP:${adapters_file}:${seed_mismatches}:${palindrome_clip_threshold}:${simple_clip_threshold}:${min_adapter_length_palindrome}:${keepbothreads}"
            } else {
                trimmomaticArgs << "ILLUMINACLIP:${adapters_file}:${seed_mismatches}:${palindrome_clip_threshold}:${simple_clip_threshold}"
            }
        }

        // Specifying the parameters for SLIDINGWINDOW
        if (window_size && required_quality) {
            trimmomaticArgs << "SLIDINGWINDOW:${window_size}:${required_quality}"
        }

        // Specifying the parameters for MAXINFO
        if (target_length && strictness) {
            trimmomaticArgs << "MAXINFO:${target_length}:${strictness}"
        }

        // Specifying the parameters for BASECOUNT
        if (bases) { 
            if (min_count) {
                if (max_count) trimmomaticArgs << "BASECOUNT:${bases}:${min_count}:${max_count}"
                else trimmomaticArgs << "BASECOUNT:${bases}:${min_count}"
            } else trimmomaticArgs << "BASECOUNT:${bases}"
        }

        // Specifying the single parameters
        if (leading) trimmomaticArgs << "LEADING:${leading}"
        if (trailing) trimmomaticArgs << "TRAILING:${trailing}"
        if (headcrop) trimmomaticArgs << "HEADCROP:${headcrop}"
        if (tailcrop) trimmomaticArgs << "TAILCROP:${tailcrop}"
        if (crop) trimmomaticArgs << "CROP:${crop}"
        if (minlen) trimmomaticArgs << "MINLEN:${minlen}"
        if (maxlen) trimmomaticArgs << "MAXLEN:${maxlen}"
        if (avgqual) trimmomaticArgs << "AVGQUAL:${avgqual}"
        
        trimmomaticArgs = trimmomaticArgs.join(" ")

        return trimmomaticArgs
    }

    static String makeBbdukArgs(qin, reads, samplerate, k, rcomp, maskmiddle, minkmerhits, minkmerfraction, mincovfraction, hammingdistance, 
    qhdist, editdistance, hammingdistance2, qhdist2, editdistance2, forbidn, ktrim, ktrimtips, maskfullycovered, mink, qtrim, trimq, minlength, 
    minlengthfraction, minavgquality, minbasequality, maxns, minconsecutivebases, trimpad, trimbyoverlap, strictoverlap, minoverlap, mininsert, 
    trimpairsevenly, forcetrimleft, forcetrimright, forcetrimright2, forcetrimmod, restrictleft, restrictright, mingc, maxgc, tossjunk) {
    
        // Add all the bbduk arguments required into the <bbdukArgs> variable
        def bbdukArgs = []

        if (qin) bbdukArgs << "qin=${qin}"
        if (reads) bbdukArgs << "reads=${reads}"
        if (samplerate) bbdukArgs << "samplerate=${samplerate}"
        if (k) bbdukArgs << "k=${k}"
        if (rcomp) bbdukArgs << "rcomp=${rcomp}"
        if (maskmiddle) bbdukArgs << "maskmiddle=${maskmiddle}"
        if (minkmerhits) bbdukArgs << "minkmerhits=${minkmerhits}"
        if (minkmerfraction) bbdukArgs << "minkmerfraction=${minkmerfraction}"
        if (mincovfraction) bbdukArgs << "mincovfraction=${mincovfraction}"
        if (hammingdistance) bbdukArgs << "hammingdistance=${hammingdistance}"
        if (qhdist) bbdukArgs << "qhdist=${qhdist}"
        if (editdistance) bbdukArgs << "editdistance=${editdistance}"
        if (hammingdistance2) bbdukArgs << "hammingdistance2=${hammingdistance2}"
        if (qhdist2) bbdukArgs << "qhdist2=${qhdist2}"
        if (editdistance2) bbdukArgs << "editdistance2=${editdistance2}"
        if (forbidn) bbdukArgs << "forbidn=${forbidn}"
        if (ktrim) bbdukArgs << "ktrim=${ktrim}"
        if (ktrimtips) bbdukArgs << "ktrimtips=${ktrimtips}"
        if (maskfullycovered) bbdukArgs << "maskfullycovered=${maskfullycovered}"
        if (mink && mink != 0) bbdukArgs << "mink=${mink}"
        if (qtrim) bbdukArgs << "qtrim=${qtrim}"
        if (trimq) bbdukArgs << "trimq=${trimq}"
        if (minlength) bbdukArgs << "minlength=${minlength}"
        if (minlengthfraction) bbdukArgs << "minlengthfraction=${minlengthfraction}"
        if (minavgquality) bbdukArgs << "minavgquality=${minavgquality}"
        if (minbasequality) bbdukArgs << "minbasequality=${minbasequality}"
        if (maxns) bbdukArgs << "maxns=${maxns}"
        if (minconsecutivebases) bbdukArgs << "minconsecutivebases=${minconsecutivebases}"
        if (trimpad) bbdukArgs << "trimpad=${trimpad}"
        if (trimbyoverlap) bbdukArgs << "trimbyoverlap=${trimbyoverlap}"
        if (strictoverlap) bbdukArgs << "strictoverlap=${strictoverlap}"
        if (minoverlap) bbdukArgs << "minoverlap=${minoverlap}"
        if (mininsert) bbdukArgs << "mininsert=${mininsert}"
        if (trimpairsevenly) bbdukArgs << "trimpairsevenly=${trimpairsevenly}"
        if (forcetrimleft) bbdukArgs << "forcetrimleft=${forcetrimleft}"
        if (forcetrimright) bbdukArgs << "forcetrimright=${forcetrimright}"
        if (forcetrimright2) bbdukArgs << "forcetrimright2=${forcetrimright2}"
        if (forcetrimmod) bbdukArgs << "forcetrimmod=${forcetrimmod}"
        if (restrictleft) bbdukArgs << "restrictleft=${restrictleft}"
        if (restrictright) bbdukArgs << "restrictright=${restrictright}"
        if (mingc) bbdukArgs << "mingc=${mingc}"
        if (maxgc) bbdukArgs << "maxgc=${maxgc}"
        if (tossjunk) bbdukArgs << "tossjunk=${tossjunk}"

        bbdukArgs = bbdukArgs.join(" ")

        return bbdukArgs
    }

    static String makeTrimGaloreArgs(paired_end, quality, quality_encoding, adapter_sequence_1, adapter_sequence_2, specific_adapters, max_length, 
    stringency, error_rate, length, maxn, trim_n, trim_1, clip_R1, clip_R2, three_prime_clip_R1, three_prime_clip_R2, nextseq_quality, hardtrim5, hardtrim3) {
        
        // If hard trimming is intended, only consider the hard trimming parameters then quit
        if (hardtrim5) {
            def hardtrimArgs = "--hardtrim5 ${hardtrim5}"
            return hardtrimArgs
        }

        if (hardtrim3) {
            def hardtrimArgs = "--hardtrim3 ${hardtrim3}"
            return hardtrimArgs
        }
        
        // Add all the trim_galore arguments required into the <trimGaloreArgs> variable
        def trimGaloreArgs = []
        if (quality_encoding) trimGaloreArgs << "--${quality_encoding}"

        // Specifying the adapters parameters
        if (adapter_sequence_1) trimGaloreArgs << "--adapter ${adapter_sequence_1}"
        if (paired_end && adapter_sequence_2) trimGaloreArgs << "--adapter2 ${adapter_sequence_2}"
        if (specific_adapters) trimGaloreArgs << "--${specific_adapters}"

        // Specify the quality
        if (nextseq_quality) trimGaloreArgs << "--nextseq ${nextseq_quality}"
        if (quality) trimGaloreArgs << "--quality ${quality}"

        // Specify the trimming parameters
        if (max_length) trimGaloreArgs << "--max_length ${max_length}"
        if (stringency) trimGaloreArgs << "--stringency ${stringency}"
        if (error_rate) trimGaloreArgs << "-e ${error_rate}"
        if (length) trimGaloreArgs << "--length ${length}"
        if (maxn) trimGaloreArgs << "--max_n ${maxn}"
        if (trim_n) trimGaloreArgs << "--trim-n"
        if (paired_end && trim_1) trimGaloreArgs << "--trim1"
        if (clip_R1) trimGaloreArgs << "--clip_R1 ${clip_R1}"
        if (paired_end && clip_R2) trimGaloreArgs << "--clip_R2 ${clip_R2}"
        if (three_prime_clip_R1) trimGaloreArgs << "--three_prime_clip_R1 ${three_prime_clip_R1}"
        if (paired_end && three_prime_clip_R2) trimGaloreArgs << "--three_prime_clip_R2 ${three_prime_clip_R2}"

        trimGaloreArgs = trimGaloreArgs.join(" ")

        return trimGaloreArgs
    }

    /* Making arguments for the alignment tools */
    static String makeStarReferenceIndexArgs(runRNGseed, genomeChrBinNbits, genomeSAindexNbases, genomeSAsparseD, genomeSuffixLengthMax, sjdbGTFfeatureExon, 
    sjdbGTFtagExonParentTranscript, sjdbGTFtagExonParentGene, sjdbGTFtagExonParentGeneName, sjdbGTFtagExonParentGeneType, sjdbOverhang, sjdbScore, sjdbInsertSave) {

        // Add all the STAR reference index arguments required into the <starReferenceIndexArgs> variable
        def starReferenceIndexArgs = []
        if (runRNGseed) starReferenceIndexArgs << "--runRNGseed ${runRNGseed}"
        if (genomeChrBinNbits) starReferenceIndexArgs << "--genomeChrBinNbits ${genomeChrBinNbits}"
        if (genomeSAindexNbases) starReferenceIndexArgs << "--genomeSAindexNbases ${genomeSAindexNbases}"
        if (genomeSAsparseD) starReferenceIndexArgs << "--genomeSAsparseD ${genomeSAsparseD}"
        if (genomeSuffixLengthMax) starReferenceIndexArgs << "--genomeSuffixLengthMax ${genomeSuffixLengthMax}"
        if (sjdbGTFfeatureExon) starReferenceIndexArgs << "--sjdbGTFfeatureExon ${sjdbGTFfeatureExon}"
        if (sjdbGTFtagExonParentTranscript) starReferenceIndexArgs << "--sjdbGTFtagExonParentTranscript ${sjdbGTFtagExonParentTranscript}"
        if (sjdbGTFtagExonParentGene) starReferenceIndexArgs << "--sjdbGTFtagExonParentGene ${sjdbGTFtagExonParentGene}"
        if (sjdbGTFtagExonParentGeneName) starReferenceIndexArgs << "--sjdbGTFtagExonParentGeneName ${sjdbGTFtagExonParentGeneName}"
        if (sjdbGTFtagExonParentGeneType) starReferenceIndexArgs << "--sjdbGTFtagExonParentGeneType ${sjdbGTFtagExonParentGeneType}"
        if (sjdbOverhang) starReferenceIndexArgs << "--sjdbOverhang ${sjdbOverhang}"
        if (sjdbScore) starReferenceIndexArgs << "--sjdbScore ${sjdbScore}"
        if (sjdbInsertSave) starReferenceIndexArgs << "--sjdbInsertSave ${sjdbInsertSave}"

        starReferenceIndexArgs = starReferenceIndexArgs.join(" ")

        return starReferenceIndexArgs  
    }

    static String makeStarArgs(paired_end, runRNGseed, readMapNumber, readMatesLengthsIn, readQualityScoreBase, clipAdapterType, clip3pNbases, 
    clip3pAdapterSeq, clip3pAdapterMMp, clip3pAfterAdapterNbases, clip5pNbases, outReadsUnmapped, outQSconversionAdd, outMultimapperOrder, outSAMmode, outSAMstrandField, 
    outSAMattributes, outSAMattrIHstart, outSAMunmapped, outSAMprimaryFlag, outSAMreadID, outSAMmapqUnique, outSAMflagOR, outSAMflagAND, outSAMmultNmax, outSAMtlen, 
    outBAMcompression, outFilterType, outFilterMultimapScoreRange, outFilterMultimapNmax, outFilterMismatchNmax, outFilterMismatchNoverLmax, outFilterMismatchNoverReadLmax, 
    outFilterScoreMin, outFilterScoreMinOverLread, outFilterMatchNmin, outFilterMatchNminOverLread, outFilterIntronMotifs, outFilterIntronStrands, outSJfilterReads, 
    outSJfilterOverhangMin, outSJfilterCountUniqueMin, outSJfilterCountTotalMin, outSJfilterDistToOtherSJmin, outSJfilterIntronMaxVsReadN, scoreGap, scoreGapNoncan, scoreGapGCAG, 
    scoreGapATAC, scoreGenomicLengthLog2scale, scoreDelOpen, scoreDelBase, scoreInsOpen, scoreInsBase, scoreStitchSJshift, seedSearchStartLmax, seedSearchStartLmaxOverLread, 
    seedSearchLmax, seedMultimapNmax, seedPerReadNmax, seedPerWindowNmax, seedNoneLociPerWindow, seedSplitMin, seedMapMin, alignIntronMin, alignIntronMax, 
    alignMatesGapMax, alignSJoverhangMin, alignSJstitchMismatchNmax, alignSJDBoverhangMin, alignSplicedMateMapLmin, alignSplicedMateMapLminOverLmate, alignWindowsPerReadNmax, 
    alignTranscriptsPerWindowNmax, alignTranscriptsPerReadNmax, alignEndsType, alignEndsProtrude, alignSoftClipAtReferenceEnds, alignInsertionFlush, peOverlapNbasesMin, 
    peOverlapMMp, winAnchorMultimapNmax, winBinNbits, winAnchorDistNbins, winFlankNbins, winReadCoverageRelativeMin, winReadCoverageBasesMin, quantMode, 
    quantTranscriptomeBAMcompression, quantTranscriptomeSAMoutput) {

        // Add all the STAR arguments required into the <starArgs> variable
        def starArgs = []
        if (runRNGseed) starArgs << "--runRNGseed ${runRNGseed}"
        if (readMapNumber) starArgs << "-readMapNumber ${readMapNumber}"
        if (readMatesLengthsIn) starArgs << "--readMatesLengthsIn ${readMatesLengthsIn}"
        if (readQualityScoreBase) starArgs << "--readQualityScoreBase ${readQualityScoreBase}"
        if (clipAdapterType) starArgs << "--clipAdapterType ${clipAdapterType}"
        if (clip3pNbases) starArgs << "--clip3pNbases ${clip3pNbases}"
        if (clip3pAdapterMMp) starArgs << "--clip3pAdapterMMp ${clip3pAdapterMMp}"
        if (clip3pAfterAdapterNbases) starArgs << "--clip3pAfterAdapterNbases ${clip3pAfterAdapterNbases}"
        if (clip5pNbases) starArgs << "--clip5pNbases ${clip5pNbases}"
        if (outReadsUnmapped) starArgs << "--outReadsUnmapped ${outReadsUnmapped}"
        if (outQSconversionAdd) starArgs << "--outQSconversionAdd ${outQSconversionAdd}"
        if (outMultimapperOrder) starArgs << "--outMultimapperOrder ${outMultimapperOrder}"
        if (outSAMmode) starArgs << "--outSAMmode ${outSAMmode}"
        if (outSAMstrandField) starArgs << "--outSAMstrandField ${outSAMstrandField}"
        if (outSAMattributes) starArgs << "--outSAMattributes ${outSAMattributes}"
        if (outSAMattrIHstart) starArgs << "--outSAMattrIHstart ${outSAMattrIHstart}"
        if (outSAMunmapped) starArgs << "--outSAMunmapped ${outSAMunmapped}"
        if (outSAMprimaryFlag) starArgs << "--outSAMprimaryFlag ${outSAMprimaryFlag}"
        if (outSAMreadID) starArgs << "--outSAMreadID ${outSAMreadID}"
        if (outSAMmapqUnique) starArgs << "--outSAMmapqUnique ${outSAMmapqUnique}"
        if (outSAMflagOR) starArgs << "--outSAMflagOR ${outSAMflagOR}"
        if (outSAMflagAND) starArgs << "--outSAMflagAND ${outSAMflagAND}"
        if (outSAMmultNmax) starArgs << "--outSAMmultNmax ${outSAMmultNmax}"
        if (outSAMtlen) starArgs << "--outSAMtlen ${outSAMtlen}"
        if (outBAMcompression) starArgs << "--outBAMcompression ${outBAMcompression}"
        if (outFilterType) starArgs << "--outFilterType ${outFilterType}"
        if (outFilterMultimapScoreRange) starArgs << "--outFilterMultimapScoreRange ${outFilterMultimapScoreRange}"
        if (outFilterMultimapNmax) starArgs << "--outFilterMultimapNmax ${outFilterMultimapNmax}"
        if (outFilterMismatchNmax) starArgs << "--outFilterMismatchNmax ${outFilterMismatchNmax}"
        if (outFilterMismatchNoverLmax) starArgs << "--outFilterMismatchNoverLmax ${outFilterMismatchNoverLmax}"
        if (outFilterMismatchNoverReadLmax) starArgs << "--outFilterMismatchNoverReadLmax ${outFilterMismatchNoverReadLmax}"
        if (outFilterScoreMin) starArgs << "--outFilterScoreMin ${outFilterScoreMin}"
        if (outFilterScoreMinOverLread) starArgs << "--outFilterScoreMinOverLread ${outFilterScoreMinOverLread}"
        if (outFilterMatchNmin) starArgs << "--outFilterMatchNmin ${outFilterMatchNmin}"
        if (outFilterMatchNminOverLread) starArgs << "--outFilterMatchNminOverLread ${outFilterMatchNminOverLread}"
        if (outFilterIntronMotifs) starArgs << "--outFilterIntronMotifs ${outFilterIntronMotifs}"
        if (outFilterIntronStrands) starArgs << "--outFilterIntronStrands ${outFilterIntronStrands}"
        if (outSJfilterReads) starArgs << "--outSJfilterReads ${outSJfilterReads}"
        if (outSJfilterOverhangMin) starArgs << "--outSJfilterOverhangMin ${outSJfilterOverhangMin}"
        if (outSJfilterCountUniqueMin) starArgs << "--outSJfilterCountUniqueMin ${outSJfilterCountUniqueMin}"
        if (outSJfilterCountTotalMin) starArgs << "--outSJfilterCountTotalMin ${outSJfilterCountTotalMin}"
        if (outSJfilterDistToOtherSJmin) starArgs << "--outSJfilterDistToOtherSJmin ${outSJfilterDistToOtherSJmin}"
        if (outSJfilterIntronMaxVsReadN) starArgs << "--outSJfilterIntronMaxVsReadN ${outSJfilterIntronMaxVsReadN}"
        if (scoreGap) starArgs << "--scoreGap ${scoreGap}"
        if (scoreGapNoncan) starArgs << "--scoreGapNoncan ${scoreGapNoncan}"
        if (scoreGapGCAG) starArgs << "--scoreGapGCAG ${scoreGapGCAG}"
        if (scoreGapATAC) starArgs << "--scoreGapATAC ${scoreGapATAC}"
        if (scoreGenomicLengthLog2scale) starArgs << "--scoreGenomicLengthLog2scale ${scoreGenomicLengthLog2scale}"
        if (scoreDelOpen) starArgs << "--scoreDelOpen ${scoreDelOpen}"
        if (scoreDelBase) starArgs << "--scoreDelBase ${scoreDelBase}"
        if (scoreInsOpen) starArgs << "--scoreInsOpen ${scoreInsOpen}"
        if (scoreInsBase) starArgs << "--scoreInsBase ${scoreInsBase}"
        if (scoreStitchSJshift) starArgs << "--scoreStitchSJshift ${scoreStitchSJshift}"
        if (seedSearchStartLmax) starArgs << "--seedSearchStartLmax ${seedSearchStartLmax}"
        if (seedSearchStartLmaxOverLread) starArgs << "--seedSearchStartLmaxOverLread ${seedSearchStartLmaxOverLread}"
        if (seedSearchLmax) starArgs << "--seedSearchLmax ${seedSearchLmax}"
        if (seedMultimapNmax) starArgs << "--seedMultimapNmax ${seedMultimapNmax}"
        if (seedPerReadNmax) starArgs << "--seedPerReadNmax ${seedPerReadNmax}"
        if (seedPerWindowNmax) starArgs << "--seedPerWindowNmax ${seedPerWindowNmax}"
        if (seedNoneLociPerWindow) starArgs << "--seedNoneLociPerWindow ${seedNoneLociPerWindow}"
        if (seedSplitMin) starArgs << "--seedSplitMin ${seedSplitMin}"
        if (seedMapMin) starArgs << "--seedMapMin ${seedMapMin}"
        if (alignIntronMin) starArgs << "--alignIntronMin ${alignIntronMin}"
        if (alignIntronMax) starArgs << "--alignIntronMax ${alignIntronMax}"
        if (alignMatesGapMax) starArgs << "--alignMatesGapMax ${alignMatesGapMax}"
        if (alignSJoverhangMin) starArgs << "--alignSJoverhangMin ${alignSJoverhangMin}"
        if (alignSJstitchMismatchNmax) starArgs << "--alignSJstitchMismatchNmax ${alignSJstitchMismatchNmax}"
        if (alignSJDBoverhangMin) starArgs << "--alignSJDBoverhangMin ${alignSJDBoverhangMin}"
        if (alignSplicedMateMapLmin) starArgs << "--alignSplicedMateMapLmin ${alignSplicedMateMapLmin}"
        if (alignSplicedMateMapLminOverLmate) starArgs << "--alignSplicedMateMapLminOverLmate ${alignSplicedMateMapLminOverLmate}"
        if (alignWindowsPerReadNmax) starArgs << "--alignWindowsPerReadNmax ${alignWindowsPerReadNmax}"
        if (alignTranscriptsPerWindowNmax) starArgs << "--alignTranscriptsPerWindowNmax ${alignTranscriptsPerWindowNmax}"
        if (alignTranscriptsPerReadNmax) starArgs << "--alignTranscriptsPerReadNmax ${alignTranscriptsPerReadNmax}"
        if (alignEndsType) starArgs << "--alignEndsType ${alignEndsType}"
        if (alignEndsProtrude) starArgs << "--alignEndsProtrude ${alignEndsProtrude}"
        if (alignSoftClipAtReferenceEnds) starArgs << "--alignSoftClipAtReferenceEnds ${alignSoftClipAtReferenceEnds}"
        if (alignInsertionFlush) starArgs << "--alignInsertionFlush ${alignInsertionFlush}"
        if (winAnchorMultimapNmax) starArgs << "--winAnchorMultimapNmax ${winAnchorMultimapNmax}"
        if (winBinNbits) starArgs << "--winBinNbits ${winBinNbits}"
        if (winAnchorDistNbins) starArgs << "--winAnchorDistNbins ${winAnchorDistNbins}"
        if (winFlankNbins) starArgs << "--winFlankNbins ${winFlankNbins}"
        if (winReadCoverageRelativeMin) starArgs << "--winReadCoverageRelativeMin ${winReadCoverageRelativeMin}"
        if (winReadCoverageBasesMin) starArgs << "--winReadCoverageBasesMin ${winReadCoverageBasesMin}"
        if (quantTranscriptomeBAMcompression) starArgs << "--quantTranscriptomeBAMcompression ${quantTranscriptomeBAMcompression}"
        if (quantTranscriptomeSAMoutput) starArgs << "--quantTranscriptomeSAMoutput ${quantTranscriptomeSAMoutput}"
        if (clip3pAdapterSeq) starArgs << "--clip3pAdapterSeq ${clip3pAdapterSeq}"
        if (quantMode) starArgs << "--quantMode ${quantMode}"

        // Add the parameters for paired-end mode
        if (paired_end) {
            if (peOverlapNbasesMin) starArgs << "--peOverlapNbasesMin ${peOverlapNbasesMin}"
            if (peOverlapMMp) starArgs << "--peOverlapMMp ${peOverlapMMp}"
        }

        starArgs = starArgs.join(" ")

        return starArgs  
    }

    static String makeHisat2ReferenceIndexArgs(large_index, noauto, bmax, bmaxdivn, dcv, nodc, noref, justref, offrate, ftabchars, localoffrate, localftabchars, seed, cutoff) {
        
        def hisat2ReferenceIndexArgs = []

        // Add the parameters specific to memory usage
        if (noauto) {
            hisat2ReferenceIndexArgs << "--noauto"
            if (bmax) hisat2ReferenceIndexArgs << "--bmax ${bmax}"
            if (bmaxdivn) hisat2ReferenceIndexArgs << "--bmaxdivn ${bmaxdivn}"
            if (dcv) hisat2ReferenceIndexArgs << "--dcv ${dcv}"
        }

        // Add all other parameters
        if (large_index) hisat2ReferenceIndexArgs << "--large-index"
        if (nodc) hisat2ReferenceIndexArgs << "--nodc"
        if (noref) hisat2ReferenceIndexArgs << "--noref"
        if (justref) hisat2ReferenceIndexArgs << "--justref"
        if (offrate) hisat2ReferenceIndexArgs << "--offrate"
        if (ftabchars) hisat2ReferenceIndexArgs << "--ftabchars ${ftabchars}"
        if (localoffrate) hisat2ReferenceIndexArgs << "--localoffrate ${localoffrate}"
        if (localftabchars) hisat2ReferenceIndexArgs << "localftabchars ${localftabchars}"
        if (seed) hisat2ReferenceIndexArgs << "--seed ${seed}"
        if (cutoff) hisat2ReferenceIndexArgs << "--cutoff ${cutoff}"

        hisat2ReferenceIndexArgs = hisat2ReferenceIndexArgs.join(" ")

        return hisat2ReferenceIndexArgs
    }

    static String makeHisat2Args(paired_end, skip, upto, trim5, trim3, phred_quality, solexa_quals, int_quals, n_ceil_func, ignore_quals, norc, nofw, mismatch_penalties, soft_clipping, 
    no_softclip, n_penalty, read_gap_penalty, reference_gap_penalty, score_min_func, pen_cansplice, pen_noncansplice, pen_canintronlen, pen_noncanintronlen, min_intronlen, 
    max_intronlen, no_temp_splicesite, no_spliced_alignment, rna_strandness, transcriptome_mapping_only, downstream_transcriptome_assembly, dta_cufflinks, avoid_pseudogene, 
    no_templatelen_adjustment, num_alignments_per_read, max_seeds, report_all_alignments, report_secondary_alignments, min_fragment_length, max_fragment_length, mate_orientations, 
    no_mixed, no_discordant, index_offrate, reorder, rng_seed, non_deterministic) {
        
        def hisat2Args = []

        // Add the parameters for paired-end mode
        if (paired_end) {
            if (min_fragment_length) hisat2Args << "--minins ${min_fragment_length}"
            if (max_fragment_length) hisat2Args << "--maxins ${max_fragment_length}"
            if (mate_orientations) hisat2Args << "--${mate_orientations}"
            if (no_mixed) hisat2Args << "--no-mixed"
            if (no_discordant) hisat2Args << "--no-discordant"          
        }

         // Add all other parameters
         // Input options
        if (skip) hisat2Args << "--skip ${skip}"
        if (upto) hisat2Args << "--upto ${upto}"
        if (trim5) hisat2Args << "--trim5 ${trim5}"
        if (trim3) hisat2Args << "--trim3 ${trim3}"
        if (phred_quality) hisat2Args << "--${phred_quality}"
        if (solexa_quals) hisat2Args << "--solexa-quals"
        if (int_quals) hisat2Args << "--int-quals"
        // Alignment options
        if (n_ceil_func) hisat2Args << "--n-ceil ${n_ceil_func}"
        if (ignore_quals) hisat2Args << "--ignore-quals"
        if (norc) hisat2Args << "--norc"
        if (nofw) hisat2Args << "--nofw"
        // Scoring options
        if (mismatch_penalties) hisat2Args << "--mp ${mismatch_penalties}"
        if (soft_clipping) hisat2Args << "--sp ${soft_clipping}"
        if (no_softclip) hisat2Args << "--no-softclip"
        if (n_penalty) hisat2Args << "--np ${n_penalty}"
        if (read_gap_penalty) hisat2Args << "--rdg ${read_gap_penalty}"
        if (reference_gap_penalty) hisat2Args << "--rfg ${reference_gap_penalty}"
        if (score_min_func) hisat2Args << "--score-min ${score_min_func}"
        // Spliced alignment options
        if (pen_cansplice) hisat2Args << "--pen-cansplice ${pen_cansplice}"
        if (pen_noncansplice) hisat2Args << "--pen-noncansplice ${pen_noncansplice}"
        if (pen_canintronlen) hisat2Args << "--pen-canintronlen ${pen_canintronlen}"
        if (pen_noncanintronlen) hisat2Args << "--pen-noncanintronlen ${pen_noncanintronlen}"
        if (min_intronlen) hisat2Args << "--min-intronlen ${min_intronlen}"
        if (max_intronlen) hisat2Args << "--max-intronlen ${max_intronlen}"
        if (no_temp_splicesite) hisat2Args << "--no-temp-splicesite"
        if (no_spliced_alignment) hisat2Args << "--no-spliced-alignment"
        if (rna_strandness) hisat2Args << "--rna-strandness ${rna_strandness}"
        if (transcriptome_mapping_only) hisat2Args << "--transcriptome-mapping-only"
        if (downstream_transcriptome_assembly) hisat2Args << "--downstream-transcriptome-assembly"
        if (dta_cufflinks) hisat2Args << "--dta-cufflinks"
        if (avoid_pseudogene) hisat2Args << "--avoid-pseudogene"
        if (no_templatelen_adjustment) hisat2Args << "--no-templatelen-adjustment"
        // Reporting options
        if (num_alignments_per_read) hisat2Args << "-k ${num_alignments_per_read}"
        if (max_seeds) hisat2Args << "--max-seeds ${max_seeds}"
        if (report_all_alignments) hisat2Args << "--all"
        if (report_secondary_alignments) hisat2Args << "--secondary"
        // Performance options
        if (index_offrate) hisat2Args << "--offrate ${index_offrate}"
        if (reorder) hisat2Args << "--reorder"
        // Other options
        if (rng_seed) hisat2Args << "--seed ${rng_seed}"
        if (non_deterministic) hisat2Args << "--non-deterministic"

        hisat2Args = hisat2Args.join(" ")

        return hisat2Args
    }

    /* Making arguments for the counting reads tools */
    static String makeHtseqCountArgs(max_reads_in_buffer, stranded, minaqual, feature_type, id_attribute, additional_attributes, mode, nonunique_mode, secondary_alignments, supplementary_alignments, 
    add_chromosome_info) {

        // Add all the HTseq arguments required into the <htseqCountArgs> variable
        def htseqCountArgs = []
        if (add_chromosome_info) htseqCountArgs << "--add-chromosome-info"
        if (max_reads_in_buffer) htseqCountArgs << "--max-reads-in-buffer=${max_reads_in_buffer}"
        if (stranded) htseqCountArgs << "--stranded=${stranded}"
        if (minaqual) htseqCountArgs << "-a ${minaqual}"
        if (feature_type) htseqCountArgs << "--type=${feature_type}"
        if (id_attribute) htseqCountArgs << "--idattr=${id_attribute}"
        if (mode) htseqCountArgs << "--mode=${mode}"
        if (nonunique_mode) htseqCountArgs << "--nonunique=${nonunique_mode}"
        if (secondary_alignments) htseqCountArgs << "--secondary-alignments=${secondary_alignments}"
        if (supplementary_alignments) htseqCountArgs << "--supplementary-alignments=${supplementary_alignments}"
        
        // Add additional attributes if specified as a list
        if (additional_attributes) {
            if (!(additional_attributes instanceof List)) throw new IllegalArgumentException("Expected a List for the additional_attributes parameter for htseq-count process, but got ${additional_attributes.getClass().name}")
            else htseqCountArgs << additional_attributes.collect { "--additional-attr=${it}" }.join(" ")
        }

        htseqCountArgs = htseqCountArgs.join(" ")

        return htseqCountArgs
    }

    static String makeFeatureCountsArgs(paired_end, requireBothEndsMapped, countChimericFragments, checkFragLength, countReadPairs, autosort, minFragLength, maxfragLength, useMetaFeatures, 
    isGTFAnnotationFile, attrType_GTF, juncCounts, isLongRead, countMultiMappingReads, allowMultiOverlap, minMQS, isStrandSpecific, featureType_GTF, byReadGroup, attrType_GTF_extra, 
    fraction, fracOverlap, fracOverlapFeature, ignoreDup, largestOverlap, minOverlap, nonOverlap, nonOverlapFeature, nonSplitOnly, primaryOnly, read2pos, readExtension3, readExtension5, 
    readShiftSize, readShiftType, splitOnly) {

        def featureCountsArgs = []

        // Add the parameters for paired-end mode
        if (paired_end) {
            if (requireBothEndsMapped) featureCountsArgs << "-B"
            if (countChimericFragments) featureCountsArgs << "-C"
            if (checkFragLength) featureCountsArgs << "-P"
            if (countReadPairs) featureCountsArgs << "--countReadPairs"
            if (autosort) featureCountsArgs << "--donotsort"
        }

        // Add the other parameters
        if (minFragLength) featureCountsArgs << "-d ${minFragLength}"
        if (maxfragLength) featureCountsArgs << "-D ${maxfragLength}"
        if (useMetaFeatures) featureCountsArgs << "-f"
        if (isGTFAnnotationFile) featureCountsArgs << "-F GTF"
        if (attrType_GTF) featureCountsArgs << "-g ${attrType_GTF}"
        if (juncCounts) featureCountsArgs << "-J"
        if (isLongRead) featureCountsArgs << "-L"
        if (countMultiMappingReads) featureCountsArgs << "-M"
        if (allowMultiOverlap) featureCountsArgs << "-O"
        if (minMQS) featureCountsArgs << "-Q ${minMQS}"
        if (isStrandSpecific) featureCountsArgs << "-s ${isStrandSpecific}"
        if (featureType_GTF) featureCountsArgs << "-t ${featureType_GTF}"
        if (byReadGroup) featureCountsArgs << "--byReadGroup"
        if (attrType_GTF_extra) featureCountsArgs << "--extraAttributes ${attrType_GTF_extra}"
        if (fraction) featureCountsArgs << "--fraction"
        if (fracOverlap) featureCountsArgs << "--fracOverlap ${fracOverlap}"
        if (fracOverlapFeature) featureCountsArgs << "--fracOverlapFeature ${fracOverlapFeature}"
        if (ignoreDup) featureCountsArgs << "--ignoreDup"
        if (largestOverlap) featureCountsArgs << "--largestOverlap"
        if (minOverlap) featureCountsArgs << "--minOverlap ${minOverlap}"
        if (nonOverlap) featureCountsArgs << "--nonOverlap ${nonOverlap}"
        if (nonOverlapFeature) featureCountsArgs << "--nonOverlapFeature ${nonOverlapFeature}"
        if (nonSplitOnly) featureCountsArgs << "--nonSplitOnly"
        if (primaryOnly) featureCountsArgs << "--primary"
        if (read2pos) featureCountsArgs << "--read2pos ${read2pos}"
        if (readExtension3) featureCountsArgs << "--readExtension3 ${readExtension3}"
        if (readExtension5) featureCountsArgs << "--readExtension5 ${readExtension5}"
        if (readShiftSize) featureCountsArgs << "--readShiftSize ${readShiftSize}"
        if (readShiftType) featureCountsArgs << "--readShiftType ${readShiftType}"
        if (splitOnly) featureCountsArgs << "--splitOnly"

        featureCountsArgs = featureCountsArgs.join(" ")

        return featureCountsArgs
    }

    static String makeSalmonMappingArgs(seqBias, gcBias, posBias, incompatPrior, meta, discardOrphansQuasi, consensusSlack, preMergeChainSubThresh, postMergeChainSubThresh, orphanChainSubThresh, 
    scoreExp, minScoreFraction, mismatchSeedSkip, disableChainingHeuristic, match_score, mismatch_score, gap_open_score, gap_extension_score, bandwidth, allowDovetail, recoverOrphans, miminBT2, 
    mimicStrictBT2, softclip, softclipOverhangs, fullLengthAlignment, hardFilter, minAlnProb, writeQualities, hitFilterPolicy, alternativeInitMode, skipQuant, dumpEq, dumpEqWeights, minAssignedFrags, 
    reduceGCMemory, biasSpeedSamp, fldMax, fldMean, fldSD, maxOccsPerHit, maxReadOcc, noLengthCorrection, noEffectiveLengthCorrection, noSingleFragProb, noFragLengthDist, noBiasLengthThreshold, 
    numBiasSamples, numAuxModelSamples, numPreAuxModelSamples, useEM, useVBOpt, rangeFactorizationBins, numGibbsSamples, noGammaDraw, numBootstraps, bootstrapReproject, thinningFactor, perTranscriptPrior, 
    perNucleotidePrior, sigDigits, vbPrior) {

        def salmonMappingArgs = []

        // Add all the salmon quant arguments required into the <salmonMappingArgs> variable
        if (seqBias) salmonMappingArgs << "--seqBias"
        if (gcBias) salmonMappingArgs << "--gcBias"
        if (posBias) salmonMappingArgs << "--posBias"
        if (incompatPrior) salmonMappingArgs << "--incompatPrior ${incompatPrior}"
        if (meta) salmonMappingArgs << "--meta"
        if (discardOrphansQuasi) salmonMappingArgs << "--discardOrphansQuasi"
        if (consensusSlack) salmonMappingArgs << "--consensusSlack ${consensusSlack}"
        if (preMergeChainSubThresh) salmonMappingArgs << "--preMergeChainSubThresh ${preMergeChainSubThresh}"
        if (postMergeChainSubThresh) salmonMappingArgs << "--postMergeChainSubThresh ${postMergeChainSubThresh}"
        if (orphanChainSubThresh) salmonMappingArgs << "--orphanChainSubThresh ${orphanChainSubThresh}"
        if (scoreExp) salmonMappingArgs << "--scoreExp ${scoreExp}"
        if (minScoreFraction) salmonMappingArgs << "--minScoreFraction ${minScoreFraction}"
        if (mismatchSeedSkip) salmonMappingArgs << "--mismatchSeedSkip ${mismatchSeedSkip}"
        if (disableChainingHeuristic) salmonMappingArgs << "--disableChainingHeuristic ${disableChainingHeuristic}"
        if (match_score) salmonMappingArgs << "--ma ${}"
        if (mismatch_score) salmonMappingArgs << "--mp ${}"
        if (gap_open_score) salmonMappingArgs << "--go ${}"
        if (gap_extension_score) salmonMappingArgs << "--ge ${}"
        if (bandwidth) salmonMappingArgs << "--bandwidth ${bandwidth}"
        if (allowDovetail) salmonMappingArgs << "--allowDovetail"
        if (recoverOrphans) salmonMappingArgs << "--recoverOrphans"
        if (miminBT2) salmonMappingArgs << "--mimicBT2"
        if (mimicStrictBT2) salmonMappingArgs << "--mimicStrictBT2"
        if (softclip) salmonMappingArgs << "--softclip"
        if (softclipOverhangs) salmonMappingArgs << "--softclipOverhangs"
        if (fullLengthAlignment) salmonMappingArgs << "--fullLengthAlignment"
        if (hardFilter) salmonMappingArgs << "--hardFilter"
        if (minAlnProb) salmonMappingArgs << "--minAlnProb ${minAlnProb}"
        if (writeQualities) salmonMappingArgs << "--writeQualities"
        if (hitFilterPolicy) salmonMappingArgs << "--hitFilterPolicy ${hitFilterPolicy}"
        if (alternativeInitMode) salmonMappingArgs << "--alternativeInitMode"
        if (skipQuant) salmonMappingArgs << "--skipQuant"
        if (dumpEq) salmonMappingArgs << "--dumpEq"
        if (dumpEqWeights) salmonMappingArgs << "--dumpEqWeights"
        if (minAssignedFrags) salmonMappingArgs << "--minAssignedFrags ${minAssignedFrags}"
        if (reduceGCMemory) salmonMappingArgs << "--reduceGCMemory"
        if (biasSpeedSamp) salmonMappingArgs << "--biasSpeedSamp ${biasSpeedSamp}"
        if (fldMax) salmonMappingArgs << "--fldMax ${fldMax}"
        if (fldMean) salmonMappingArgs << "--fldMean ${fldMean}"
        if (fldSD) salmonMappingArgs << "--fldSD ${fldSD}"
        if (maxOccsPerHit) salmonMappingArgs << "--maxOccsPerHit ${maxOccsPerHit}"
        if (maxReadOcc) salmonMappingArgs << "--maxReadOcc ${maxReadOcc}"
        if (noLengthCorrection) salmonMappingArgs << "--noLengthCorrection"
        if (noEffectiveLengthCorrection) salmonMappingArgs << "--noEffectiveLengthCorrection"
        if (noSingleFragProb) salmonMappingArgs << "--noSingleFragProb"
        if (noFragLengthDist) salmonMappingArgs << "--noFragLengthDist"
        if (noBiasLengthThreshold) salmonMappingArgs << "--noBiasLengthThreshold" 
        if (numBiasSamples) salmonMappingArgs << "--numBiasSamples ${numBiasSamples}"
        if (numAuxModelSamples) salmonMappingArgs << "--numAuxModelSamples ${numAuxModelSamples}"
        if (numPreAuxModelSamples) salmonMappingArgs << "--numPreAuxModelSamples ${numPreAuxModelSamples}"
        if (useEM) salmonMappingArgs << "--useEM"
        if (useVBOpt) salmonMappingArgs << "--useVBOpt"
        if (rangeFactorizationBins) salmonMappingArgs << "--rangeFactorizationBins ${rangeFactorizationBins}"
        if (numGibbsSamples) salmonMappingArgs << "--numGibbsSamples ${numGibbsSamples}"
        if (noGammaDraw) salmonMappingArgs << "--noGammaDraw"
        if (numBootstraps) salmonMappingArgs << "--numBootstraps ${numBootstraps}"
        if (bootstrapReproject) salmonMappingArgs << "--bootstrapReproject"
        if (thinningFactor) salmonMappingArgs << "--thinningFactor ${thinningFactor}"
        if (perTranscriptPrior) salmonMappingArgs << "--perTranscriptPrior"
        if (perNucleotidePrior) salmonMappingArgs << "--perNucleotidePrior"
        if (sigDigits) salmonMappingArgs << "--sigDigits ${sigDigits}"
        if (vbPrior) salmonMappingArgs << "--vbPrior ${vbPrior}"

        salmonMappingArgs = salmonMappingArgs.join(" ")

        return salmonMappingArgs
    }
 

}

