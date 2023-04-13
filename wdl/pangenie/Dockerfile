FROM condaforge/mambaforge:latest as conda

ADD . /pangenie
WORKDIR /pangenie

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    pkg-config \
    samtools \
    tabix \
    bcftools && \
    apt-get -y clean  && \
    apt-get -y autoclean  && \
    apt-get -y autoremove

RUN conda env create -f environment.yml
RUN echo "source activate pangenie" > ~/.bashrc

SHELL ["/bin/bash", "--login", "-c"]
RUN which jellyfish
RUN mkdir /pangenie/build && cd /pangenie/build && cmake .. && make

CMD ["/bin/bash", "--login", "-c"]
