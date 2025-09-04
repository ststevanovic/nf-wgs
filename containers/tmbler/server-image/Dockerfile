# docker image aguida/tmbler
#
# usage
#   > docker run -it aguida/tmbler
#
#   This  will launch an interactive R console with the package already loaded
#
# to build the container
#   $ docker build -t aguida/tmbler:latest .
#
#   This is the bash command to execute to make the build of the Dockerfile
#
# to run it interactively in bash
#   $ docker run -it aguida/tmbler:latest bash
#
#   This is the terminal command to enter inside the Docker container interactively
#   in bash instead of R as described in the docker file.
# ------------------------------------

FROM aguida/tmbler:latest

COPY . /usr/local/src/tmbler
WORKDIR /usr/local/src/tmbler

RUN apt-get update

# Install Git lfs to download large file
RUN apt-get install -y curl
RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash
RUN apt-get install -y git-lfs

# This is a dependency for Roxygen2 to work as expeceted
RUN apt-get install -y qpdf
# this is a dependency for some of the plotting libraries
RUN sudo apt install -y libmagick++-dev

RUN sudo apt-get install -y libbz2-dev
RUN sudo apt-get install -y liblzma-dev

# Install dependencies
RUN R -e "install.packages('devtools')"
# RUN R -e "devtools::install('TMBleR', dependencies=TRUE)"

WORKDIR /scripts
COPY /scripts/tmbler.R .
