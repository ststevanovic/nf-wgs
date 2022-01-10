#!/bin/bash

# BUILDS
FASTQC_IMAGE_NAME="fastqc:v0.11.9"
FASTQC_CONTAINER_DIR="fastqc"

SAMTOOLS_IMAGE_NAME="samtools:v0.11.9"
SAMTOOLS_CONTAINER_DIR="samtools"

MULTIQC_IMAGE_NAME="multiqc:latest"
MULTIQC_CONTAINER_DIR="multiqc"

cd $FASTQC_CONTAINER_DIR
docker build -t $FASTQC_IMAGE_NAME .
cd ..

cd $GATK_CONTAINER_DIR
docker build -t $GATK_IMAGE_NAME .
cd ..

cd $SAMTOOLS_CONTAINER_DIR
docker build -t $SAMTOOLS_IMAGE_NAME .
cd ..

cd $MULTIQC_CONTAINER_DIR
docker build -t $MULTIQC_IMAGE_NAME .
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





