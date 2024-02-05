FROM condaforge/mambaforge:22.11.1-4 as conda

ADD . /kage-lite
WORKDIR /kage-lite

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends \
    build-essential \
    samtools \
    tabix \
    bcftools \
    zlib1g-dev && \
    apt-get -y clean  && \
    apt-get -y autoclean  && \
    apt-get -y autoremove

RUN pip install numpy==1.23.5 pandas==2.1.1 truvari==4.1.0 scikit-learn==1.3.1 matplotlib==3.8.0 pytest==8.0.0
RUN pip install -e shared_memory_wrapper/ -e obgraph/ -e graph_kmer_index/ -e kmer_mapper/ -e kage/ -e lite_utils/

CMD ["bash"]