FROM condaforge/mambaforge:latest as conda

ADD . /kage-lite
WORKDIR /kage-lite

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends \
    build-essential \
    samtools \
    tabix \
    bcftools \
    zlib1g && \
    apt-get -y clean  && \
    apt-get -y autoclean  && \
    apt-get -y autoremove

RUN pip install pandas truvari
RUN pip install -e shared_memory_wrapper/ -e obgraph/ -e graph_kmer_index/ -e kmer_mapper/ -e kage/ -e lite_utils/

CMD ["bash"]