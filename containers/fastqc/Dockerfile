################# BASE IMAGE ######################
FROM biocontainers/biocontainers:v1.1.0_cv2

################## METADATA ######################
LABEL \
    base_image="biocontainers:v1.1.0_cv2" \
    version="8" \
    software="fastqc" \
    software.version="0.11.9" \
    about.summary="A quality control tool for high throughput sequence data." \
    about.home="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/" \
    about.documentation="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/" \
    about.license_file="https://www.gnu.org/copyleft/gpl.html" \
    about.license="SPDX:GPL-3.0" \
    extra.identifiers.biotools="fastqc" \
    about.tags="General" \
    extra.binaries="fastqc" \
    about.maintainer="Alvin Chen <ychen@aperiomics.com>" 

################## INSTALLATION ######################
USER root

ENV ZIP=fastqc_v0.11.9.zip
ENV URL=https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
ENV FOLDER=FastQC
ENV DST=/opt/fastqc
ENV ULOCAL=/usr/local/bin

RUN mkdir -p $DST
RUN apt-get install build-essential
RUN apt-get update && apt-get install -y openjdk-8-jre-headless 
RUN wget $URL/$ZIP -O $DST/$ZIP --no-check-certificate
RUN unzip $DST/$ZIP -d $DST 
RUN rm $DST/$ZIP 
RUN mv $DST/$FOLDER/* $DST 
RUN rmdir $DST/$FOLDER 
RUN chmod +x $DST/fastqc 
RUN ln -s $DST/fastqc $ULOCAL/fastqc

WORKDIR /data/
# fastqc binary is available in path, as a wrapper to java command
# graphical interface needs X11, not available in container, only 
# command-line should be used (non-interactive)