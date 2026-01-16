class OrganizeArguments {

    /* Making arguments for the trimming tools */
    static String makeTrimmomaticArgs(base_quality_encoding, adapters_file, seed_mismatches, palindrome_clip_threshold, simple_clip_threshold, 
    min_adapter_length_palindrome, keepbothreads, window_size, required_quality, target_length, strictness, bases, 
    min_count, max_count, leading, trailing, headcrop, tailcrop, crop, minlen, maxlen, avgqual) {

        // Add all the trimmomatic arguments required into the <trimmomaticArgs> variable
        def trimmomaticArgs = ["-${base_quality_encoding}"]

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
        def bbdukArgs = [
            "qin=${qin}", "reads=${reads}", "samplerate=${samplerate}", "k=${k}", "rcomp=${rcomp}", "maskmiddle=${maskmiddle}", 
            "minkmerhits=${minkmerhits}", "minkmerfraction=${minkmerfraction}", "mincovfraction=${mincovfraction}", "hammingdistance=${hammingdistance}", 
            "qhdist=${qhdist}", "editdistance=${editdistance}", "hammingdistance2=${hammingdistance2}", "qhdist2=${qhdist2}", 
            "editdistance2=${editdistance2}", "forbidn=${forbidn}", "ktrim=${ktrim}", "ktrimtips=${ktrimtips}", "maskfullycovered=${maskfullycovered}", 
            "mink=${mink}", "qtrim=${qtrim}", "trimq=${trimq}", "minlength=${minlength}", "minlengthfraction=${minlengthfraction}", 
            "minavgquality=${minavgquality}", "minbasequality=${minbasequality}", "maxns=${maxns}", "minconsecutivebases=${minconsecutivebases}", 
            "trimpad=${trimpad}", "trimbyoverlap=${trimbyoverlap}", "strictoverlap=${strictoverlap}", "minoverlap=${minoverlap}", 
            "mininsert=${mininsert}", "trimpairsevenly=${trimpairsevenly}", "forcetrimleft=${forcetrimleft}", "forcetrimright=${forcetrimright}", 
            "forcetrimright2=${forcetrimright2}", "forcetrimmod=${forcetrimmod}", "restrictleft=${restrictleft}", "restrictright=${restrictright}", 
            "mingc=${mingc}", "maxgc=${maxgc}", "tossjunk=${tossjunk}"
        ]

        if (mink == 0) bbdukArgs.remove("mink=${mink}")
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
        else trimGaloreArgs << "--quality ${quality}"

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
        def starReferenceIndexArgs = ["--runRNGseed ${runRNGseed}", "--genomeChrBinNbits ${genomeChrBinNbits}", "--genomeSAindexNbases ${genomeSAindexNbases}", "--genomeSAsparseD ${genomeSAsparseD}", 
            "--genomeSuffixLengthMax ${genomeSuffixLengthMax}", "--sjdbGTFfeatureExon ${sjdbGTFfeatureExon}", "--sjdbGTFtagExonParentTranscript ${sjdbGTFtagExonParentTranscript}", 
            "--sjdbGTFtagExonParentGene ${sjdbGTFtagExonParentGene}", "--sjdbGTFtagExonParentGeneName ${sjdbGTFtagExonParentGeneName}", "--sjdbGTFtagExonParentGeneType ${sjdbGTFtagExonParentGeneType}", 
            "--sjdbOverhang ${sjdbOverhang}", "--sjdbScore ${sjdbScore}", "--sjdbInsertSave ${sjdbInsertSave}"]

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
        def starArgs = ["--runRNGseed ${runRNGseed}", "--readMapNumber ${readMapNumber}", "--readMatesLengthsIn ${readMatesLengthsIn}", "--readQualityScoreBase ${readQualityScoreBase}", "--clipAdapterType ${clipAdapterType}", "--clip3pNbases ${clip3pNbases}", 
            "--clip3pAdapterMMp ${clip3pAdapterMMp}", "--clip3pAfterAdapterNbases ${clip3pAfterAdapterNbases}", "--clip5pNbases ${clip5pNbases}", "--outReadsUnmapped ${outReadsUnmapped}", "--outQSconversionAdd ${outQSconversionAdd}", 
            "--outMultimapperOrder ${outMultimapperOrder}", "--outSAMmode ${outSAMmode}", "--outSAMstrandField ${outSAMstrandField}", "--outSAMattributes ${outSAMattributes}", "--outSAMattrIHstart ${outSAMattrIHstart}", 
            "--outSAMunmapped ${outSAMunmapped}", "--outSAMprimaryFlag ${outSAMprimaryFlag}", "--outSAMreadID ${outSAMreadID}", "--outSAMmapqUnique ${outSAMmapqUnique}", "--outSAMflagOR ${outSAMflagOR}", 
            "--outSAMflagAND ${outSAMflagAND}", "--outSAMmultNmax ${outSAMmultNmax}", "--outSAMtlen ${outSAMtlen}", "--outBAMcompression ${outBAMcompression}", "--outFilterType ${outFilterType}", 
            "--outFilterMultimapScoreRange ${outFilterMultimapScoreRange}", "--outFilterMultimapNmax ${outFilterMultimapNmax}", "--outFilterMismatchNmax ${outFilterMismatchNmax}", "--outFilterMismatchNoverLmax ${outFilterMismatchNoverLmax}", 
            "--outFilterMismatchNoverReadLmax ${outFilterMismatchNoverReadLmax}", "--outFilterScoreMin ${outFilterScoreMin}", "--outFilterScoreMinOverLread ${outFilterScoreMinOverLread}", "--outFilterMatchNmin ${outFilterMatchNmin}", 
            "--outFilterMatchNminOverLread ${outFilterMatchNminOverLread}", "--outFilterIntronMotifs ${outFilterIntronMotifs}", "--outFilterIntronStrands ${outFilterIntronStrands}", "--outSJfilterReads ${outSJfilterReads}", 
            "--outSJfilterOverhangMin ${outSJfilterOverhangMin}", "--outSJfilterCountUniqueMin ${outSJfilterCountUniqueMin}", "--outSJfilterCountTotalMin ${outSJfilterCountTotalMin}", "--outSJfilterDistToOtherSJmin ${outSJfilterDistToOtherSJmin}", 
            "--outSJfilterIntronMaxVsReadN ${outSJfilterIntronMaxVsReadN}", "--scoreGap ${scoreGap}", "--scoreGapNoncan ${scoreGapNoncan}", "--scoreGapGCAG ${scoreGapGCAG}", "--scoreGapATAC ${scoreGapATAC}", 
            "--scoreGenomicLengthLog2scale ${scoreGenomicLengthLog2scale}", "--scoreDelOpen ${scoreDelOpen}", "--scoreDelBase ${scoreDelBase}", "--scoreInsOpen ${scoreInsOpen}", "--scoreInsBase ${scoreInsBase}", 
            "--scoreStitchSJshift ${scoreStitchSJshift}", "--seedSearchStartLmax ${seedSearchStartLmax}", "--seedSearchStartLmaxOverLread ${seedSearchStartLmaxOverLread}", "--seedSearchLmax ${seedSearchLmax}", 
            "--seedMultimapNmax ${seedMultimapNmax}", "--seedPerReadNmax ${seedPerReadNmax}", "--seedPerWindowNmax ${seedPerWindowNmax}", "--seedNoneLociPerWindow ${seedNoneLociPerWindow}", "--seedSplitMin ${seedSplitMin}", 
            "--seedMapMin ${seedMapMin}", "--alignIntronMin ${alignIntronMin}", "--alignIntronMax ${alignIntronMax}", "--alignMatesGapMax ${alignMatesGapMax}", "--alignSJoverhangMin ${alignSJoverhangMin}", 
            "--alignSJstitchMismatchNmax ${alignSJstitchMismatchNmax}", "--alignSJDBoverhangMin ${alignSJDBoverhangMin}", "--alignSplicedMateMapLmin ${alignSplicedMateMapLmin}", "--alignSplicedMateMapLminOverLmate ${alignSplicedMateMapLminOverLmate}", 
            "--alignWindowsPerReadNmax ${alignWindowsPerReadNmax}", "--alignTranscriptsPerWindowNmax ${alignTranscriptsPerWindowNmax}", "--alignTranscriptsPerReadNmax ${alignTranscriptsPerReadNmax}", "--alignEndsType ${alignEndsType}", 
            "--alignEndsProtrude ${alignEndsProtrude}", "--alignSoftClipAtReferenceEnds ${alignSoftClipAtReferenceEnds}", "--alignInsertionFlush ${alignInsertionFlush}", "--winAnchorMultimapNmax ${winAnchorMultimapNmax}", 
            "--winBinNbits ${winBinNbits}", "--winAnchorDistNbins ${winAnchorDistNbins}", "--winFlankNbins ${winFlankNbins}", "--winReadCoverageRelativeMin ${winReadCoverageRelativeMin}", "--winReadCoverageBasesMin ${winReadCoverageBasesMin}", 
            "--quantTranscriptomeBAMcompression ${quantTranscriptomeBAMcompression}", "--quantTranscriptomeSAMoutput ${quantTranscriptomeSAMoutput}"]
        
        // Add the optional parameters
        if (clip3pAdapterSeq) starArgs << "--clip3pAdapterSeq ${clip3pAdapterSeq}"
        if (quantMode) starArgs << "--quantMode ${quantMode}"

        // Add the parameters for paired-end mode
        if (paired_end) {
            starArgs << "--peOverlapNbasesMin ${peOverlapNbasesMin}"
            starArgs << "--peOverlapMMp ${peOverlapMMp}"
        }

        starArgs = starArgs.join(" ")

        return starArgs  
    }

    /* Making arguments for the counting reads tools */
    static String makeHtseqCountArgs(max_reads_in_buffer, stranded, minaqual, feature_type, id_attribute, additional_attributes, mode, nonunique_mode, secondary_alignments, supplementary_alignments, 
    add_chromosome_info) {

        // Add all the HTseq arguments required into the <htseqCountArgs> variable
        def htseqCountArgs = ["--max-reads-in-buffer=${max_reads_in_buffer}", "--stranded=${stranded}", "-a ${minaqual}", "--type=${feature_type}", "--idattr=${id_attribute}", "--mode=${mode}", 
            "--nonunique=${nonunique_mode}", "--secondary-alignments=${secondary_alignments}", "--supplementary-alignments=${supplementary_alignments}"]
        
        // Add additional attributes if specified as a list
        if (additional_attributes) {
            if (!(additional_attributes instanceof List)) throw new IllegalArgumentException("Expected a List for the additional_attributes parameter for htseq-count process, but got ${additional_attributes.getClass().name}")
            else htseqCountArgs << additional_attributes.collect { "--additional-attr=${it}" }.join(" ")
        }

        // Check other parameters
        if (add_chromosome_info) htseqCountArgs << "--add-chromosome-info"

        htseqCountArgs = htseqCountArgs.join(" ")

        return htseqCountArgs

    }



}
