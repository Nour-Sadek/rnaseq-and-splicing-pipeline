class OrganizeArguments {

    static String makeTrimmomaticArgs(base_quality_encoding, adapters_file, seed_mismatches, palindrome_clip_threshold, simple_clip_threshold, 
    min_adapter_length_palindrome, keepbothreads, window_size, required_quality, target_length, strictness, bases, 
    min_count, max_count, leading, trailing, headcrop, tailcrop, crop, minlen, maxlen, avgqual) {

        // Add all the trimmomatic arguments required into the <trimmomaticArgs> variable
        def trimmomaticArgs = ["-${base_quality_encoding}"]

        // Specifying the parameters for ILLUMINACLIP
        if (seed_mismatches != "none" && palindrome_clip_threshold != "none" && simple_clip_threshold != "none") {
            if (min_adapter_length_palindrome != "none" && keepbothreads != "none") {
                trimmomaticArgs << "ILLUMINACLIP:${adapters_file}:${seed_mismatches}:${palindrome_clip_threshold}:${simple_clip_threshold}:${min_adapter_length_palindrome}:${keepbothreads}"
            } else {
                trimmomaticArgs << "ILLUMINACLIP:${adapters_file}:${seed_mismatches}:${palindrome_clip_threshold}:${simple_clip_threshold}"
            }
        }

        // Specifying the parameters for SLIDINGWINDOW
        if (window_size != "none" && required_quality != "none") {
            trimmomaticArgs << "SLIDINGWINDOW:${window_size}:${required_quality}"
        }

        // Specifying the parameters for MAXINFO
        if (target_length != "none" && strictness != "none") {
            trimmomaticArgs << "MAXINFO:${target_length}:${strictness}"
        }

        // Specifying the parameters for BASECOUNT
        if (bases != "none") { 
            if (min_count != "none") {
                if (max_count != "none") trimmomaticArgs << "BASECOUNT:${bases}:${min_count}:${max_count}"
                else trimmomaticArgs << "BASECOUNT:${bases}:${min_count}"
            } else trimmomaticArgs << "BASECOUNT:${bases}"
        }

        // Specifying the single parameters
        if (leading != "none") trimmomaticArgs << "LEADING:${leading}"
        if (trailing != "none") trimmomaticArgs << "TRAILING:${trailing}"
        if (headcrop != "none") trimmomaticArgs << "HEADCROP:${headcrop}"
        if (tailcrop != "none") trimmomaticArgs << "TAILCROP:${tailcrop}"
        if (crop != "none") trimmomaticArgs << "CROP:${crop}"
        if (minlen != "none") trimmomaticArgs << "MINLEN:${minlen}"
        if (maxlen != "none") trimmomaticArgs << "MAXLEN:${maxlen}"
        if (avgqual != "none") trimmomaticArgs << "AVGQUAL:${avgqual}"
        
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
        if (hardtrim5 != 'none') {
            def hardtrimArgs = "--hardtrim5 ${hardtrim5}"
            return hardtrimArgs
        }

        if (hardtrim3 != 'none') {
            def hardtrimArgs = "--hardtrim3 ${hardtrim3}"
            return hardtrimArgs
        }
        
        // Add all the trim_galore arguments required into the <trimGaloreArgs> variable
        def trimGaloreArgs = []
        if (quality_encoding != 'none') trimGaloreArgs << "--${quality_encoding}"

        // Specifying the adapters parameters
        if (adapter_sequence_1 != 'none') trimGaloreArgs << "--adapter ${adapter_sequence_1}"
        if (paired_end && adapter_sequence_2 != 'none') trimGaloreArgs << "--adapter2 ${adapter_sequence_2}"
        if (specific_adapters != 'none') trimGaloreArgs << "--${specific_adapters}"

        // Specify the quality
        if (nextseq_quality != 'none') trimGaloreArgs << "--nextseq ${nextseq_quality}"
        else trimGaloreArgs << "--quality ${quality}"

        // Specify the trimming parameters
        if (max_length != 'none') trimGaloreArgs << "--max_length ${max_length}"
        if (stringency != 'none') trimGaloreArgs << "--stringency ${stringency}"
        if (error_rate != 'none') trimGaloreArgs << "-e ${error_rate}"
        if (length != 'none') trimGaloreArgs << "--length ${length}"
        if (maxn != 'none') trimGaloreArgs << "--max_n ${maxn}"
        if (trim_n) trimGaloreArgs << "--trim-n"
        if (paired_end && trim_1) trimGaloreArgs << "--trim1"
        if (clip_R1 != 'none') trimGaloreArgs << "--clip_R1 ${clip_R1}"
        if (paired_end && clip_R2 != 'none') trimGaloreArgs << "--clip_R2 ${clip_R2}"
        if (three_prime_clip_R1 != 'none') trimGaloreArgs << "--three_prime_clip_R1 ${three_prime_clip_R1}"
        if (paired_end && three_prime_clip_R2 != 'none') trimGaloreArgs << "--three_prime_clip_R2 ${three_prime_clip_R2}"

        trimGaloreArgs = trimGaloreArgs.join(" ")

        return trimGaloreArgs
    }
    
}
