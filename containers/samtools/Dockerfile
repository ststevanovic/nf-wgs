FROM ubuntu:latest 

ENV \
    DST=/opt/samtools \
    VER=1.14

LABEL \
    version=$VER \
    software="samtools" \
    software.version="0.11.9" \
    about.summary="A a suite of programs for interacting with high-throughput sequencing data." \
    about.home="http://www.htslib.org/" \
    about.documentation="http://www.htslib.org/doc/" \
    about.license_file="https://www.gnu.org/copyleft/gpl.html" \
    about.license="CC BY-NC-ND 4.0" \
    about.tags="General"


RUN apt-get update && apt install curl -y

RUN apt-get install build-essential libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev -y
RUN curl -sSL https://github.com/samtools/samtools/releases/download/$VER/samtools-$VER.tar.bz2 | tar xj
RUN cd samtools-$VER \ 
    && ./configure --prefix=$DST \
    && make \
    && make install 
ENV PATH="$DST/bin:$PATH"
# RUN export PATH="$DST/bin:${PATH}"
# RUN cd /root && rm -r samtools-$VER 