# Description of the parameters for each of the tools represented in this pipeline

The information is essentially copied from the original documents for the tools, but it felt appropriate to have a summary of all of the parameters of the tools here, in one place.

## Trimmmomatic

Please refer to the original github page for more information: https://github.com/usadellab/Trimmomatic

- `base_quality_encoding`: Specifies the quality score encoding. Options are: 'phred64' and 'phred33'. phred33 is the standard for modern Illumina data.
- `adapters_file`: used in ILLUMINACLIP command. Specifies which fasta file to use that has all the adapters, PCR sequences, etc. Trimmomatic already has multiple fasta files that contain Illumina adapter and other technical sequences. Options are: 'NexteraPE-PE.fa', 'TruSeq2-PE.fa', 'TruSeq2-SE.fa', 'TruSeq3-PE-2-GGGGG.fa', "TruSeq3-PE-2.fa", 'TruSeq3-PE.fa', 'TruSeq3-SE.fa'
- `seed_mismatches`: used in ILLUMINACLIP command. Specifies the maximum mismatch count which will still allow a full match to be performed.
- `palindrome_clip_threshold`: used in ILLUMINACLIP command. Specifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment.
- `simple_clip_threshold`: used in ILLUMINACLIP command. Specifies how accurate the match between any adapter etc. sequence must be against a read.
- `min_adapter_length_palindrome`: used in ILLUMINACLIP command. (optional, if the other three ILLUMINACLIP command parameters were specified) Specifies the minimum adapter length in palindrome mode [default = 8].
- `keepbothreads`: used in ILLUMINACLIP command. (optional boolean) specifies if both reads should be kept in palindrome mode even when redundant information is found (small inserts)[default = False]. Note: min_adapter_length_palindrome needs to be set manually to be able to activate keepbothreads.
- `window_size`: used in SLIDINGWINDOW command. Specifies the number of bases to average across.
- `required_quality`: used in SLIDINGWINDOW command. Specifies the average quality required.
- `leading`: used in LEADING command. Specifies the minimum quality required to keep a base.
- `trailing`: used in TRAILING command. Specifies the minimum quality required to keep a base.
- `headcrop`: used in HEADCROP command. Specifies the number of bases to remove from the start of the read.
- `tailcrop`: used in TAILCROP command. Specifies the number of bases to remove from the end of the read.
- `crop`: used in CROP command. Specifies the number of bases to keep, from the start of the read.
- `target_length`: used in MAXINFO command. Specifies the ideal length for a read. Scoring system favors reads trimmed to a length near this value.
- `strictness`: used in MAXINFO command. Specifies the trade-off between length and quality. A higher strictness value places more weight on base quality and less on the length score. Range between 0.0 and 1.0.
- `minlen`: used in MINLEN command. Specifies the minimum length of reads to be kept.
- `maxlen`: used in MAXLEN command. Specifies the maximum length of reads to be kept.
- `avgqual`: used in AVGQUAL command. Specifies the minimum Phred quality score of a read to be kept.
- `bases`: used in BASECOUNT command. Specifies the string of characters to be counted (e.g., N, GC).
- `min_count`: used in BASECOUNT command. (optional, if bases was specified) Specifies the minimum number of times the bases must appear for the read to be kept.
- `max_count`: used in BASECOUNT command. (optional, if bases and min_count parameters of the BASECOUNT command were specified) Specifies the maximum number of times the bases are allowed to appear.

## BBDUK

Please refer to the official website for more information: https://bbmap.org/tools/bbduk

Input parameters:
- `qin`: Input quality offset: 33 (Sanger), 64, or auto.
- `reads`: If positive, quit after processing X reads or pairs.
- `samplerate`: Default: 1. Set lower to only process a fraction of input reads.

Processing parameters:
- `k`: Kmer length used for finding contaminants. Contaminants shorter than k will not be found. k must be at least 1. For k>31, BBDuk uses emulated long kmers by requiring consecutive 31-mer matches.
- `rcomp`: Look for reverse-complements of kmers in addition to forward kmers. Set rcomp=f to match only forward orientation if exact directional matching is required.
- `maskmiddle`: (mm) Default: t. Treat the middle base of a kmer as a wildcard, to increase sensitivity in the presence of errors. This may also be set to a number, e.g. mm=3, to mask that many bp. Default mm=t corresponds to mm=1 for odd-length kmers and mm=2 for even-length kmers, while mm=f equals mm=0. Disabled when using mink or kbig>k.
- `minkmerhits`: (mkh) Default: 1. Reads need at least this many matching kmers to be considered as matching the reference.
- `minkmerfraction`: (mkf) Default: 0.0. A read needs at least this fraction of its total kmers to hit a ref, in order to be considered a match. If this and minkmerhits are set, the greater is used.
- `mincovfraction`: (mcf) Default: 0.0. A read needs at least this fraction of its total bases to be covered by ref kmers to be considered a match. If specified, mcf overrides mkh and mkf.
- `hammingdistance`: (hdist) Deafult: 0. Maximum Hamming distance for ref kmers (subs only). Memory use scales as (3×K)^hdist: E. coli with k=31, hdist=1 uses ~15GB vs 140MB at hdist=0.
- `qhdist`: Default: 0. Hamming distance for query kmers; impacts speed, not memory. Alternative to hdist for memory-constrained systems - mutates read kmers instead of reference kmers.
- `editdistance`: (edist) Default: 0. Maximum edit distance from ref kmers (subs and indels). Memory use scales as (8×K)^edist. Use qhdist instead of edist/hdist if insufficient memory available.
- `hammingdistance2`: (hdist2) Deafult: 0. Sets hdist for short kmers, when using mink.
- `qhdist2`: Default: 0. Sets qhdist for short kmers, when using mink.
- `editdistance2`: (edist2) Default: 0. Sets edist for short kmers, when using mink.
- `forbidn`: (fn) Default: f. Forbids matching of read kmers containing N. By default, these will match a reference 'A' if hdist>0 or edist>0, to increase sensitivity.

Trimming/filtering/masking parameters:
- `ktrim`: Default: f. Trim reads to remove bases matching reference kmers, plus all bases to the left or right. Values: f (don't trim), r (trim to the right), l (trim to the left)
- `ktrimtips`: Default: 0. Set this to a positive number to perform ktrim on both ends, examining only the outermost X bases.
- `maskfullycovered`: (mfc) Default: f. Only mask bases that are fully covered by kmers.
- `mink`: Default: 0. Look for shorter kmers at read tips down to this length, when k-trimming or masking. 0 means disabled. Essential for adapter trimming when adapter remnants are shorter than k. Enabling this will disable maskmiddle.
- `qtrim`: Default: f. Trim read ends to remove bases with quality below trimq. Performed AFTER looking for kmers. Uses Phred algorithm for more accurate trimming than naive approaches. Values: rl (trim both ends), f (neither end), r (right end only), l (left end only), w (sliding window).
- `trimq`: Default: 6. Regions with average quality BELOW this will be trimmed, if qtrim is set to something other than f. Can be a floating-point number like 7.3.
- `minlength`: (ml) Default: 10. Reads shorter than this after trimming will be discarded. Pairs will be discarded if both are shorter.
- `minlengthfraction`: (mlf) Default: 0. Reads shorter than this fraction of original length after trimming will be discarded.
- `minavgquality`: (maq) Reads with average quality (after trimming) below this will be discarded.
- `minbasequality`: (mbq) Default: 0. Reads with any base below this quality (after trimming) will be discarded.
- `maxns`: Default: -1. If non-negative, reads with more Ns than this (after trimming) will be discarded.
- `minconsecutivebases`: (mcb) Default: 0. Discard reads without at least this many consecutive called bases.
- `trimpad`: (tp) Default: 0. Trim this much extra around matching kmers.
- `trimbyoverlap`: (tbo) Default: f. Trim adapters based on where paired reads overlap. Uses BBMerge internally and does not require known adapter sequences.
- `strictoverlap`: Default: t. Adjust sensitivity for trimbyoverlap mode.
- `minoverlap`: Default: 14. Require this many bases of overlap for detection.
- `mininsert`: Default: 40. Require insert size of at least this for overlap. Should be reduced to 16 for small RNA sequencing.
- `trimpairsevenly`: (tpe) Default: f. When kmer right-trimming, trim both reads to the minimum length of either. Useful when adapter kmer detected in only one read of a pair.
- `forcetrimleft`: (ftl) Default: 0. If positive, trim bases to the left of this position (exclusive, 0-based).
- `forcetrimright`: (ftr) Default: 0. If positive, trim bases to the right of this position (exclusive, 0-based).
- `forcetrimright2`: (ftr2) Default: 0.  If positive, trim this many bases on the right end.
- `forcetrimmod`: (ftm) Default: 0. If positive, right-trim length to be equal to zero, modulo this number. Use ftm=5 to remove inaccurate extra bases from 151bp reads (converts to 150bp), common in Illumina runs.
- `restrictleft`: Default: 0. If positive, only look for kmer matches in the leftmost X bases.
- `restrictright`: Default: 0. If positive, only look for kmer matches in the rightmost X bases.
- `mingc`: Default: 0. Discard reads with GC content below this.
- `maxgc`: Default: 1. Discard reads with GC content above this.
- `tossjunk`: Default: f. Discard reads with invalid characters as bases.

## Trim Galore

Please refer to the official user guide on their GitHub page for more information: https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md

- `quality`: Default: 20. Trim low-quality ends from reads in addition to adapter removal in a single pass.
- `quality_encoding`: Instructs Cutadapt to use ASCII+33 quality scores as Phred scores for phred33 and ASCII+64 for phred64 for quality trimming.
- `adapter_sequence_1`: Adapter sequence to be trimmed. If not specified explicitly, Trim Galore will try to auto-detect whether the Illumina universal, Nextera transposase or Illumina small RNA adapter sequence was used. If no adapter can be detected within the first 1 million sequences of the first file specified Trim Galore defaults to --illumina. A single base may also be given as e.g. -a A{10}, to be expanded to -a AAAAAAAAAA.
- `adapter_sequence_2`: Optional adapter sequence to be trimmed off read 2 of paired-end files. This option requires --paired to be specified as well.
- `specific_adapters`: Options are `illumina` for the adapter sequence to be trimmed be the first 13bp of the illumina universal adapter `AGATCGGAAGAGC`, `stranded_illumina` for the first 13bp of the Illumina stranded mRNA or Total RNA adapter `ACTGTCTCTTATA`, `nextera` for the first 12bp of the Nextera adapter `CTGTCTCTTATA`, `small_rna` for the first 12bp of the Illumina Small RNA 3' Adapter `TGGAATTCTCGG`, and `bgiseq` for BGISEQ/DNBSEQ/MGISEQ instead the default auto-detection (uses sequences AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA for Read 1 (BGI/MGI forward), and AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG for Read 2 (BGI/MGI reverse)).
- `max_length`: Discard reads that are longer than bp after trimming. This is only advised for smallRNA sequencing to remove non-small RNA sequences.
- `stringency`: Default: 1. Overlap with adapter sequence required to trim a sequence.
- `error_rate`: Default: 0.1. Maximum allowed error rate (no. of errors divided by the length of the matching region).
- `length`: Default: 20. Discard reads that became shorter than length INT because of either quality or adapter trimming. A value of 0 effectively disables this behaviour.
- `maxn`: The total number of Ns (as integer) a read may contain before it will be removed altogether.
- `trim_n`: Removes Ns from either side of the read.
- `trim_1`: Trims 1 bp off every read from its 3' end. Only used for paired-end reads
- `clip_R1`: Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads). This may be useful if the qualities were very poor, or if there is some sort of unwanted bias at the 5' end.
- `clip_R2`: Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only). This may be useful if the qualities were very poor, or if there is some sort of unwanted bias at the 5' end.
- `three_prime_clip_R1`: Instructs Trim Galore to remove <int> bp from the 3' end of read 1 (or single-end reads) AFTER adapter/quality trimming has been performed. This may remove some unwanted bias from the 3' end that is not directly related to adapter sequence or basecall quality.
- `three_prime_clip_R2`: Instructs Trim Galore to remove bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed. This may remove some unwanted bias from the 3' end that is not directly related to adapter sequence or basecall quality.
- `nextseq_quality`: This enables the option --nextseq-trim=3'CUTOFF within Cutadapt, which will set a quality cutoff (that is normally given with -quality instead), but qualities of G bases are ignored. This trimming is in common for the NextSeq- and NovaSeq-platforms, where basecalls without any signal are called as high-quality G bases. This is mutually exclusive with -quality (if this is specified, then -quality determined by the `quality` parameter)
- `hardtrim5`: Instead of performing adapter-/quality trimming, this option will simply hard-trim sequences to bp from the 3'-end. Once hard-trimming of files is complete, Trim Galore will exit.
- `harftrim3`: Instead of performing adapter-/quality trimming, this option will simply hard-trim sequences to bp from the 5'-end. Once hard-trimming of files is complete, Trim Galore will exit.
