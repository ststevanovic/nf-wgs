################## BASE IMAGE ######################
FROM ubuntu:20.04.3 

################## METADATA ######################
LABEL  \
    base_image="ubuntu:20.04.3" \
    software="bwa" \
    software.version="0.7.17" \
    about.summary="Burrows Weeler Aligner is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome." \
    about.home="http://bio-bwa.sourceforge.net/" \
    about.documentation="http://bio-bwa.sourceforge.net/bwa.shtml" \
    about.license_file="http://bio-bwa.sourceforge.net/bwa.shtml#13" \
    about.license="GPLv3" \
    about.tags="Genomics" 

################## INSTALLATION ######################
ENV VERSION=0.7.17
RUN apt-get update && apt-get install make && apt -y install gcc
RUN apt-get install libz-dev

RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

RUN wget \
    https://sourceforge.net/projects/bio-bwa/files/bwa-$VERSION.tar.bz2 \
    && tar xvf bwa-$VERSION.tar.bz2 

WORKDIR /bwa-$VERSION/
RUN make 

ENV PATH=$PATH:$PWD/bwa-$VERSION

WORKDIR / 

CMD ["bwa"]