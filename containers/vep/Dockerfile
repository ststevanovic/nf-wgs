FROM nfcore/base:1.12.1
LABEL \
    author="Maxime Garcia" \
    description="VEP image for use in nf-core/sarek" \
    maintainer="maxime.garcia@scilifelab.se"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-sarek-vep-2.7.1/bin:$PATH

# Setup default ARG variables
ARG GENOME=GRCh38
ARG SPECIES=homo_sapiens
ARG VEP_VERSION=99

# Download Genome
RUN vep_install \
    -a c \
    -c .vep \
    -s ${SPECIES} \
    -y ${GENOME} \
    --CACHE_VERSION ${VEP_VERSION} \
    --CONVERT \
    --NO_BIOPERL --NO_HTSLIB --NO_TEST --NO_UPDATE

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-sarek-vep-2.7.1 > nf-core-sarek-vep-2.7.1.yml
