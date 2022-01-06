FROM debian:stretch-slim

RUN apt-get update && apt install curl

RUN curl -sSL https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/salmon-1.5.2_linux_x86_64.tar.gz | tar xz \
    && mv /salmon-*/bin/* /usr/bin/ \
    && mv /salmon-*/lib/* /usr/lib/

ENV PATH=$PATH:/usr/games/