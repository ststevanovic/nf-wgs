################## BASE IMAGE ######################
FROM ubuntu:latest

################## METADATA ######################
LABEL  \
    base_image="ubuntu:latest" \
    version="2" \
    software="bwa-mem2" \
    software.version="2.2.1" \
    about.summary="Burrow-Wheeler Aligner for pairwise alignment between DNA sequences" \
    about.home="https://github.com/bwa-mem2/bwa-mem2" \
    about.documentation="http://bio-bwa.sourceforge.net/" \
    about.license_file="https://anaconda.org/bioconda/bwa-mem2" \
    about.license="MIT" \
    about.tags="Genomics" \
    extra.identifiers.biotools="bwa-mem2"

################## INSTALLATION ######################
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN apt-get update

RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 

RUN conda install -c bioconda bwa-mem2 -y

CMD ["bwa-mem2"]