# Pull and transform docker to apptainer containers
# It follows this format: apptainer pull --name <docker_container_name.img> docker://<docker_container_name> > /dev/null

cd $NXF_APPTAINER_CACHEDIR
apptainer pull --name community.wave.seqera.io-library-fastqc-0.12.1--af7a5314d5015c29.img docker://community.wave.seqera.io/library/fastqc:0.12.1--af7a5314d5015c29 > /dev/null
apptainer pull --name biocontainers-bowtie2-v2.4.1_cv1.img docker://biocontainers/bowtie2:v2.4.1_cv1 > /dev/null
apptainer pull --name community.wave.seqera.io-library-trimmomatic-0.40--0c25090769939729.img docker://community.wave.seqera.io/library/trimmomatic:0.40--0c25090769939729 > /dev/null
apptainer pull --name community.wave.seqera.io-library-trim-galore-0.6.10--1bf8ca4e1967cd18.img docker://community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18 > /dev/null
apptainer pull --name community.wave.seqera.io-library-bbmap-39.33--60639c9e1473b7a8.img docker://community.wave.seqera.io/library/bbmap:39.33--60639c9e1473b7a8 > /dev/null
apptainer pull --name community.wave.seqera.io-library-star-2.7.11b--822039d47adf19a7.img docker://community.wave.seqera.io/library/star:2.7.11b--822039d47adf19a7 > /dev/null
apptainer pull --name community.wave.seqera.io-library-hisat2-2.2.1--df34d2bb25ac6de5.img docker://community.wave.seqera.io/library/hisat2:2.2.1--df34d2bb25ac6de5 > /dev/null
apptainer pull --name community.wave.seqera.io-library-samtools-1.22.1--eccb42ff8fb55509.img docker://community.wave.seqera.io/library/samtools:1.22.1--eccb42ff8fb55509 > /dev/null
apptainer pull --name community.wave.seqera.io-library-htseq-2.0.9--4a65a9021e1142a5.img docker://community.wave.seqera.io/library/htseq:2.0.9--4a65a9021e1142a5 > /dev/null
apptainer pull --name community.wave.seqera.io-library-subread-2.1.1--0ac4d7e46cd0c5d7.img docker://community.wave.seqera.io/library/subread:2.1.1--0ac4d7e46cd0c5d7 > /dev/null
apptainer pull --name community.wave.seqera.io-library-salmon-1.10.3--fcd0755dd8abb423.img docker://community.wave.seqera.io/library/salmon:1.10.3--fcd0755dd8abb423 > /dev/null
apptainer pull --name community.wave.seqera.io-library-kallisto-0.51.1--b63691b6841c7a52.img docker://community.wave.seqera.io/library/kallisto:0.51.1--b63691b6841c7a52 > /dev/null
apptainer pull --name mcfonsecalab-majiq-2.5.1.img docker://mcfonsecalab/majiq:2.5.1 > /dev/null
apptainer pull --name mcfonsecalab-rmats.img docker://mcfonsecalab/rmats > /dev/null
apptainer pull --name naotokubota-suppa-2.3.img docker://naotokubota/suppa:2.3 > /dev/null
apptainer pull --name naotokubota-whippet-1.6.1.img docker://naotokubota/whippet:1.6.1 > /dev/null
apptainer pull --name rocker-r-base-4.3.2.img docker://rocker/r-base:4.3.2 > /dev/null
