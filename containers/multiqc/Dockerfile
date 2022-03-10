################## BASE IMAGE ######################
FROM ubuntu:20.04

################## METADATA ######################
LABEL  \
    base_image="ubuntu" \
    version="20.04.3" \
    software="multiqc" \
    software.version="1.11" \
    about.summary="Aggregate results from bioinformatics analyses across many samples into a single report." \
    about.home="https://github.com/ewels/MultiQC" \
    about.license="GPL-3.0" 

################## INSTALLATION ######################
RUN apt-get update --assume-yes
RUN apt install python3-pip --assume-yes
RUN pip install multiqc

CMD ["multiqc"]