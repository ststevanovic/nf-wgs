#!/bin/bash

# BUILDS
FASTQC_IMAGE_NAME="fastqc:v0.11.9"
FASTQC_CONTAINER_DIR="fastqc"

BWASAM_IMAGE_NAME="bwa-sam:latest"
BWASAM_CONTAINER_DIR="bwa-sam"

MULTIQC_IMAGE_NAME="multiqc:v1.11"
MULTIQC_CONTAINER_DIR="multiqc"

MANTA_IMAGE_NAME="manta:v1.6.0"
MANTA_CONTAINER="manta"

SNPEFF_IMAGE_NAME="snpeff:v4.3"
SNPEFF_CONTAINER="snpeff"

cd $FASTQC_CONTAINER_DIR
docker build -t $FASTQC_IMAGE_NAME .
cd ..

cd $BWA_CONTAINER_DIR
docker build -t $BWA_IMAGE_NAME .
cd ..

cd $BWASAM_CONTAINER_DIR
docker build -t $BWASAM_IMAGE_NAME .
cd ..

cd $MULTIQC_CONTAINER_DIR
docker build -t $MULTIQC_IMAGE_NAME .
cd ..

cd $MANTA_CONTAINER_DIR
docker build -t $MANTA_IMAGE_NAME .
cd ..

cd $SNPEFF_CONTAINER_DIR
docker build -t $SNPEFF_IMAGE_NAME .
cd ..

# HUB
GATK_IMAGE_NAME="gatk:latest"
GATK_CONTAINER="gatk"

QUALIMAP_IMAGE_NAME="qualimap:latest"
QUALIMAP_CONTAINER="qualimap"

cd $GATK_CONTAINER_DIR
docker pull broadinstitute/gatk
docker image tag broadinstitute/gatk:latest $GATK_IMAGE_NAME
docker rmi broadinstitute/gatk
cd ..

cd $QUALIMAP_CONTAINER
docker pull pegi3s/qualimap
docker image tag pegi3s/qualimap:latest $QUALIMAP_IMAGE_NAME
docker rmi pegi3s/qualimap
cd ..




