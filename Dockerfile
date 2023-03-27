FROM condaforge/mambaforge:latest as conda

ADD . /kage-lite
WORKDIR /kage-lite

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends \
    build-essential && \
    apt-get -y clean  && \
    apt-get -y autoclean  && \
    apt-get -y autoremove

RUN pip install -e shared_memory_wrapper/ obgraph/ graph_kmer_index/ kage/ lite_utils/

CMD ["bash"]