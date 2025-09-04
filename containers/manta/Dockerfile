FROM ubuntu:20.04
RUN apt-get -y update
RUN apt-get install -y wget \
        python

RUN mkdir manta
RUN wget --no-check-certificate https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2
RUN tar xvjf manta-1.6.0.centos6_x86_64.tar.bz2 -C manta --strip-components=1

ENV PATH=/manta/bin:$PATH