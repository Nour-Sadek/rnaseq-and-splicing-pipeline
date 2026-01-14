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
    
}
