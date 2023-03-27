FROM condaforge/mambaforge:latest as conda

ADD . /kage-lite
WORKDIR /kage-lite

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends \
    build-essential \
    samtools \
    tabix \
    bcftools && \
    apt-get -y clean  && \
    apt-get -y autoclean  && \
    apt-get -y autoremove

RUN pip install shared_memory_wrapper/ obgraph/ graph_kmer_index/ -e kage/ lite_utils/

CMD ["bash"]