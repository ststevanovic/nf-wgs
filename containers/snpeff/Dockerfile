FROM nfcore/base:1.12.1

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/snpeff-4.3/bin:$PATH

# Setup default ARG variables
ARG GENOME=GRCh38
ARG SNPEFF_CACHE_VERSION=86

# Download Genome
RUN snpEff download ${GENOME}.${SNPEFF_CACHE_VERSION}

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name snpeff-4.3 > snpeff-4.3.yml
