################## BASE IMAGE ######################
FROM ubuntu:20.04


ENV \
    BWA_VERSION=0.7.17 \
    SAMTOOLS_DST=/opt/samtools \
    SAMTOOLS_VERSION=1.14 

################## METADATA ######################
LABEL  \
    base_image="ubuntu:20.04" \
    software="bwa and samtools" \
    software.version_bwa=$BWA_VERSION \
    software.version_samtools=$SAMTOOLS_VERSION \
    about.summary="Combination of BWA and SAMTOOLS in docker file." \
    about.documentation_bwa="http://bio-bwa.sourceforge.net/bwa.shtml" \
    about.license_file_bwa="http://bio-bwa.sourceforge.net/bwa.shtml#13" \
    about.documentation_bwa="http://bio-bwa.sourceforge.net/bwa.shtml" \
    about.documentation_samtools="http://www.htslib.org/doc/" \
    about.license_file_bwa="http://bio-bwa.sourceforge.net/bwa.shtml#13" \
    about.license_file_samtools="https://www.gnu.org/copyleft/gpl.html" \
    about.license_bwa="GPLv3" \
    about.license_samtools="CC BY-NC-ND 4.0" \
    about.tags="Genomics"  

################## INSTALLATION ######################

RUN apt-get update && apt-get install make && apt -y install gcc
RUN apt-get install libz-dev build-essential libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev -y

RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

RUN wget \
    https://sourceforge.net/projects/bio-bwa/files/bwa-$BWA_VERSION.tar.bz2 \
    && tar xvf bwa-$BWA_VERSION.tar.bz2 

WORKDIR /bwa-$BWA_VERSION/
RUN make 
ENV PATH=$PATH:$PWD/bwa-$BWA_VERSION

WORKDIR / 

RUN wget -qO- \
    https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VERSION/samtools-$SAMTOOLS_VERSION.tar.bz2 | tar xfj -

WORKDIR /samtools-$SAMTOOLS_VERSION
RUN ./configure --prefix=$SAMTOOLS_DST \
    && make \
    && make install 
ENV PATH="$SAMTOOLS_DST/bin:$PATH"

WORKDIR / 

CMD ["bwa","samtools"]