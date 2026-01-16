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
- `qin`: Default: auto. Input quality offset: 33 (Sanger), 64, or auto.
- `reads`: Default: -1. If positive, quit after processing X reads or pairs.
- `samplerate`: Default: 1. Set lower to only process a fraction of input reads.

Processing parameters:
- `k`: Default: 27. Kmer length used for finding contaminants. Contaminants shorter than k will not be found. k must be at least 1. For k>31, BBDuk uses emulated long kmers by requiring consecutive 31-mer matches.
- `rcomp`: Default: t. Look for reverse-complements of kmers in addition to forward kmers. Set rcomp=f to match only forward orientation if exact directional matching is required.
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
- `specific_adapters`: Default: illumina. Options are `illumina` for the adapter sequence to be trimmed be the first 13bp of the illumina universal adapter `AGATCGGAAGAGC`, `stranded_illumina` for the first 13bp of the Illumina stranded mRNA or Total RNA adapter `ACTGTCTCTTATA`, `nextera` for the first 12bp of the Nextera adapter `CTGTCTCTTATA`, `small_rna` for the first 12bp of the Illumina Small RNA 3' Adapter `TGGAATTCTCGG`, and `bgiseq` for BGISEQ/DNBSEQ/MGISEQ instead the default auto-detection (uses sequences AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA for Read 1 (BGI/MGI forward), and AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG for Read 2 (BGI/MGI reverse)).
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

## Star

Please refer to the official user guide on their GitHub page for more information: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

Run parameters:
- `runRNGseed`: Default: 777. Random number generator seed.

Genome Indexing Parameters - only used with –runMode genomeGenerate:
- `genomeChrBinNbits`: Default: 18. int: =log2(chrBin), where chrBin is the size of the bins for genome storage: each chromosome will occupy an integer number of bins.
- `genomeSAindexNbases`: Default: 12. int: length (bases) of the SA pre-indexing string. Typically between 10 and 15. Longer strings will use much more memory, but allow faster searches. For small genomes, the parameter –genomeSAindexNbases must be scaled down to min(14, log2(GenomeLength)/2 - 1).
- `genomeSAsparseD`: Default: 1. Suffux array sparsity, i.e. distance between indices: use bigger numbers to decrease needed RAM at the cost of mapping speed reduction
- `genomeSuffixLengthMax`: Default: -1. Maximum length of the suffixes, has to be longer than read length. -1 = infinite.

Splice Junctions Database:
- `sjdbGTFfeatureExon`: Default: exon. Feature type in GTF file to be used as exons for building transcripts.
- `sjdbGTFtagExonParentTranscript`: Default: transcript_id. GTF attribute name for parent transcript ID.
- `sjdbGTFtagExonParentGene`: Default: gene_id. GTF attribute name for parent gene ID.
- `sjdbGTFtagExonParentGeneName`: Default: gene_name. GTF attribute name for parent gene name.
- `sjdbGTFtagExonParentGeneType`: Default: gene_type gene_biotype. GTF attribute name for parent gene type.
- `sjdbOverhang`: Default: 100.  Length of the donor/acceptor sequence on each side of the junctions, ideally = (mate length - 1).
- `sjdbScore`: Default: 2. Extra alignment score for alignments that cross database junctions.
- `sjdbInsertSave`: Default: Basic. Which files to save when sjdb junctions are inserted on the fly at the mapping step.
- `readMapNumber`: Default: -1. Number of reads to map from the beginning of the file.
- `readMatesLengthsIn`: Default: NotEqual. Equal/NotEqual - lengths of names,sequences,qualities for both mates are the same / not the same. NotEqual is safe in all situations.
- `readQualityScoreBase`: Default: 33. number to be subtracted from the ASCII code to get Phred quality score

Read Clipping:
- `clipAdapterType`: Default: Hamming. Adapter clipping type. Options are: Hamming, CellRanger4, or None.
- `clip3pNbases`: Default: 0.  Number(s) of bases to clip from 3p of each mate. If one value is given, it will be assumed the same for both mates.
- `clip3pAdapterSeq`: Default: null. Adapter sequences to clip from 3p of each mate. If one value is given, it will be assumed the same for both mates.
- `clip3pAdapterMMp`: Default: 0.1. Max proportion of mismatches for 3p adapter clipping for each mate. If one value is given, it will be assumed the same for both mates.
- `clip3pAfterAdapterNbases`: Default: 0. Number of bases to clip from 3p of each mate after the adapter clipping. If one value is given, it will be assumed the same for both mates.
- `clip5pNbases`: Default: 0. Number(s) of bases to clip from 5p of each mate. If one value is given, it will be assumed the same for both mates.

Output - General:
- `outReadsUnmapped`: Default: None. Output of unmapped and partially mapped (i.e. mapped only one mate of a paired end read) reads in separate file(s). Options are: None or Fastx.
- `outQSconversionAdd`: Default: 0. Add this number to the quality score (e.g. to convert from Illumina to Sanger, use -31).
- `outMultimapperOrder`: Default: Old_2.4. Order of multimapping alignments in the output files. Options are: Old_2.4 or Random.

Output - SAM and BAM:
- `outSAMmode`: Default: Full. Mode of SAM output. Options are: None, Full, or NoQS.
- `outSAMstrandField`: Default: None. Cufflinks-like strand field flag. Options are: None or intronMotif.
- `outSAMattributes`: Default: Standard. A string of desired SAM attributes, in the order desired for the output SAM. Tags can be listed in any combination/order. Many options are available (check the manual for the full list).
- `outSAMattrIHstart`: Default: 1. Start value for the IH attribute. 0 may be required by some downstream software, such as Cufflinks or StringTie.
- `outSAMunmapped`: Default: None. Output of unmapped reads in the SAM format.
- `outSAMprimaryFlag`: Default: OneBestScore. Which alignments are considered primary - all others will be marked with 0x100 bit in the FLAG
- `outSAMreadID`: Default: Standard. Read ID record type. Options are: Standard or Number.
- `outSAMmapqUnique`: Default: 255. The MAPQ value for unique mappers. Options are values 0 to 255.
- `outSAMflagOR`: Default: 0. sam FLAG will be bitwise OR’d with this value, i.e. FLAG=FLAG — outSAMflagOR. This is applied after all flags have been set by STAR, and after outSAMflagAND. Can be used to set specific bits that are not set otherwise. Options are values 0 to 65535.
- `outSAMflagAND`: Default: 65535. sam FLAG will be bitwise AND’d with this value, i.e. FLAG=FLAG & outSAMflagOR. This is applied after all flags have been set by STAR, but before outSAMflagOR. Can be used to unset specific bits that are not set otherwise. Options are values 0 to 65535.
- `outSAMmultNmax`: Default: -1. max number of multiple alignments for a read that will be output to the SAM/BAM files. Note that if this value is not equal to -1, the top scoring alignment will be output first.
- `outSAMtlen`: Default: 1. Calculation method for the TLEN field in the SAM/BAM files. Options are: 1 or 2.
- `outBAMcompression`: Default: 1. BAM compression level, -1=default compression (6?), 0=no compression, 10=maximum compression. Options are values -1 to 10.

Output Filtering:
- `outFilterType`: Default: Normal. Tyoe of filtering. Options are: Normal or BySJout.
- `outFilterMultimapScoreRange`: Default: 1. The score range below the maximum score for multimapping alignments.
- `outFilterMultimapNmax`: Default: 10. maximum number of loci the read is allowed to map to. Alignments (all of them) will be output only if the read maps to no more loci than this value. Otherwise no alignments will be output, and the read will be counted as "mapped to too many loci" in the Log.final.out.
- `outFilterMismatchNmax`: Default: 10. Alignment will be output only if it has no more mismatches than this value.
- `outFilterMismatchNoverLmax`: Default: 0.3. Alignment will be output only if its ratio of mismatches to mapped length is less than or equal to this value.
- `outFilterMismatchNoverReadLmax`: Default: 1.0. Alignment will be output only if its ratio of mismatches to read length is less than or equal to this value.
- `outFilterScoreMin`: Default: 0. Alignment will be output only if its score is higher than or equal to this value.
- `outFilterScoreMinOverLread`: Default: 0.66. Same as outFilterScoreMin, but normalized to read length (sum of mates' lengths for paired-end reads).
- `outFilterMatchNmin`: Default: 0.  Alignment will be output only if the number of matched bases is higher than or equal to this value.
- `outFilterMatchNminOverLread`: Default: 0.66. Same as outFilterMatchNmin, but normalized to the read length (sum of mates' lengths for paired-end reads).
- `outFilterIntronMotifs`: Default: None. Filter alignment using their motifs. Options are: None, RemoveNoncanonical, or RemoveNoncanonicalUnannotated.
- `outFilterIntronStrands`: Default: RemoveInconsistentStrands. Filter strands. Options are: RemoveInconsistentStrands or None.

Output Filtering - Splice Junctions:
- `outSJfilterReads`: Default: All. Which reads to consider for collapsed splice junctions output. Options are: All or Unique.
- `outSJfilterOverhangMin`: Default: 30 12 12 12. 4 integers: minimum overhang length for splice junctions on both sides for: (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif. -1 means no output for that motif. Does not apply to annotatoed junctions.
- `outSJfilterCountUniqueMin`: Default: 3 1 1 1. 4 integers: minimum uniquely mapping read count per junction for: (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif. -1 means no output for that motif. Junctions are output if one of outSJfilterCountUniqueMin OR outSJfilterCountTotalMin conditions are satisfied. Does not apply to annotatoed junctions.
- `outSJfilterCountTotalMin`: Default: 3 1 1 1. 4 integers: minimum total (multi-mapping+unique) read count per junction for: (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif. -1 means no output for that
motif. Junctions are output if one of outSJfilterCountUniqueMin OR outSJfilterCountTotalMin conditions are satisfied. Does not apply to annotatoed junctions.
- `outSJfilterDistToOtherSJmin`: Default: 10 0 5 10. 4 integers>=0: minimum allowed distance to other junctions’ donor/acceptor. Does not apply to annotatoed junctions.
- `outSJfilterIntronMaxVsReadN`: Default: 50000 100000 200000. N integers>=0: maximum gap allowed for junctions supported by 1,2,3,,,N reads i.e. by default junctions supported by 1 read can have gaps <=50000b, by 2 reads: <=100000b, by 3 reads: <=200000. by >=4 reads any gap <=alignIntronMax. Does not apply to annotatoed junctions.

Scoring:
- `scoreGap`: Default: 0. Splice junction penalty (independent on intron motif).
- `scoreGapNoncan`: Default: -8. Non-canonical junction penalty (in addition to scoreGap).
- `scoreGapGCAG`: Default: -4. GC/AG and CT/GC junction penalty (in addition to scoreGap).
- `scoreGapATAC`: Default: -8. AT/AC and GT/AT junction penalty (in addition to scoreGap).
- `scoreGenomicLengthLog2scale`: Default: -0.25. Extra score logarithmically scaled with genomic length of the alignment: scoreGenomicLengthLog2scale*log2(genomicLength).
- `scoreDelOpen`: Default: -2. Deletion open penalty.
- `scoreDelBase`: Default: -2.  Deletion extension penalty per base (in addition to scoreDelOpen).
- `scoreInsOpen`: Default: -2. Insertion open penalty.
- `scoreInsBase`: Default: -2. Insertion extension penalty per base (in addition to scoreInsOpen).
- `scoreStitchSJshift`: Default: 1. Maximum score reduction while searching for SJ boundaries in the stitching step.

Alignments and Seeding:
- `seedSearchStartLmax`: Default: 50. Defines the search start point through the read - the read is split into pieces no longer than this value.
- `seedSearchStartLmaxOverLread`: Default: 1.0. seedSearchStartLmax normalized to read length (sum of mates' lengths for paired-end reads).
- `seedSearchLmax`: Default: 0. Defines the maximum length of the seeds, if =0 seed length is not limited.
- `seedMultimapNmax`: Default: 10000. Only pieces that map fewer than this value are utilized in the stitching procedure.
- `seedPerReadNmax`: Default: 1000. Max number of seeds per read.
- `seedPerWindowNmax`: Default: 50. Max number of seeds per window.
- `seedNoneLociPerWindow`: Default: 10. Max number of one seed loci per window.
- `seedSplitMin`: Default: 12. Min length of the seed sequences split by Ns or mate gap.
- `seedMapMin`: Default: 5. Min length of seeds to be mapped.
- `alignIntronMin`: Default: 21. Minimum intron size, genomic gap is considered intron if its length>=alignIntronMin, otherwise it is considered Deletion.
- `alignIntronMax`: Default: 0. Maximum intron size, if 0, max intron size will be determined by (2ˆwinBinNbits)*winAnchorDistNbins.
- `alignMatesGapMax`: Default: 0. Maximum gap between two mates, if 0, max intron gap will be determined by (2ˆwinBinNbits)*winAnchorDistNbins.
- `alignSJoverhangMin`: Default: 5. Minimum overhang (i.e. block size) for spliced alignments.
- `alignSJstitchMismatchNmax`: Default: 0 -1 0 0. 4*int>=0: maximum number of mismatches for stitching of the splice junctions (-1: no limit). (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif.
- `alignSJDBoverhangMin`: Default: 3. Minimum overhang (i.e. block size) for annotated (sjdb) spliced.
- `alignSplicedMateMapLmin`: Default: 0. Minimum mapped length for a read mate that is spliced.
- `alignSplicedMateMapLminOverLmate`: Default: 0.66. alignSplicedMateMapLmin normalized to mate length.
- `alignWindowsPerReadNmax`: Default: 10000. Max number of windows per read.
- `alignTranscriptsPerWindowNmax`: Default: 100. Max number of transcripts per window.
- `alignTranscriptsPerReadNmax`: Default: 10000. Max number of different alignments per read to consider.
- `alignEndsType`: Default: Local. Type of read ends alignment. Options are: Local, EndToEnd, Extend5pOfRead1, or Extend5pOfReads12.
- `alignEndsProtrude`: Default: 0 ConcordantPair. Allow protrusion of alignment ends, i.e. start (end) of the +strand mate downstream of the start (end) of the -strand mate. 1st word: maximum number of protrusion bases allowed. 2nd word: Options are either ConcordantPair or DiscordantPair.
- `alignSoftClipAtReferenceEnds`: Default: Yes. Allow the soft-clipping of the alignments past the end of the chromosomes. Options are: Yes or No.
- `alignInsertionFlush`: Default: None. How to flush ambiguous insertion positions. Options are: None or Right.

Paired-End reads:
- `peOverlapNbasesMin`: Default: 0. Minimum number of overlapping bases to trigger mates merging and realignment. Specify >0 value to switch on the "merginf of overlapping mates" algorithm.
- `peOverlapMMp`: Default: 0.01. Maximum proportion of mismatched bases in the overlap area. Options are real values between 0 and 1.

Windows, Anchors, Binning:
- `winAnchorMultimapNmax`: Default: 50. Max number of loci anchors are allowed to map to.
- `winBinNbits`: Default: 16. =log2(winBin), where winBin is the size of the bin for the windows/clustering, each window will occupy an integer number of bins.
- `winAnchorDistNbins`: Default: 9. Max number of bins between two anchors that allows aggregation of anchors into one window.
- `winFlankNbins`: Default: 4. log2(winFlank), where win Flank is the size of the left and right flanking regions for each window.
- `winReadCoverageRelativeMin`: Default: 0.5. Minimum relative coverage of the read sequence by the seeds in a window, for STARlong algorithm only.
- `winReadCoverageBasesMin`: Default: 0. Minimum number of bases covered by the seeds in a window, for STARlong algorithm only.
- `quantMode`: Default: null. Types of quantification requested. If you want to use it, options are either TranscriptomeSAM or GeneCounts.
- `quantTranscriptomeBAMcompression`: Default: 1. Transcriptome BAM compression level. Options are values between -2 to 10.
- `quantTranscriptomeSAMoutput`: Default: BanSingleEnd_BanIndels_ExtendSoftclip. Alignment filtering for TranscriptomeSAM output. Options are: BanSingleEnd_BanIndels_ExtendSoftclip, BanSingleEnd_BanIndels_ExtendSoftclip, BanSingleEnd, or BanSingleEnd_ExtendSoftclip.

## HTSeq

Please refer to the official website for more information: https://htseq.readthedocs.io/en/latest/htseqcount.html#htseqcount

- `max_reads_in_buffer`: Default: 30000000. When alignment_file is paired end sorted by position, allow only so many reads to stay in memory until the mates are found (raising this number will use more memory). Has no effect for single end or paired end sorted by name.
- `stranded`: Default: 'no'. Whether the data is from a strand-specific assay. For stranded=no, a read is considered overlapping with a feature regardless of whether it is mapped to the same or the opposite strand as the feature. For stranded=yes and single-end reads, the read has to be mapped to the same strand as the feature. For paired-end reads, the first read has to be on the same strand and the second read on the opposite strand. For stranded=reverse, these rules are reversed. Options are: yes, no, or reverse.
- `minaqual`: Default: 10. Skip all reads with MAPQ alignment quality lower than the given minimum value. MAPQ is the 5th column of a SAM/BAM file and its usage depends on the software used to map the reads.
- `feature_type`: Default: 'exon'. Feature type (3rd column in GTF file) to be used, all features of other type are ignored (default, suitable for RNA-Seq analysis using an Ensembl GTF file: exon).
- `id_attribute`: Default: 'gene_id'. GTF attribute to be used as feature ID. Several GTF lines with the same feature ID will be considered as parts of the same feature. The feature ID is used to identity the counts in the output table. The default, suitable for RNA-Seq analysis using an Ensembl GTF file, is gene_id.
- `additional_attributes`: Default: null. Additional feature attributes, which will be printed as an additional column after the primary attribute column but before the counts column(s). The default is none, a suitable value to get gene names using an Ensembl GTF file is gene_name. To specify these attributes, present them as a list of string, e.g. ['gene_name', 'exon_number'].
- `mode`: Default: 'union'. Mode to handle reads overlapping more than one feature. Options are: union, intersection-strict, or intersection-nonempty.
- `nonunique_mode`: Default: 'none'. Mode to handle reads that align to or are assigned to more than one feature in the overlap mode of choice. Options are: none, all, fraction, or random.
- `secondary_alignments`: Default: 'score'. Mode to handle secondary alignments (SAM flag 0x100). Options are: score or ignore.
- `supplementary_alignments`: Default: 'score'. Mode to handle supplementary/chimeric alignments (SAM flag 0x800). Options are: score or ignore.
- `add_chromosome_info`: Default: false. Store information about the chromosome of each feature as an additional attribute (e.g. column in the TSV output file).
