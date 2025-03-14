FROM condaforge/mambaforge:22.11.1-4 as conda

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends \
    build-essential \
    tabix \
    bcftools \
    autoconf \
    zlib1g-dev \
    libncurses-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4 \
    libcrypto++6 && \
    apt-get -y clean && \
    apt-get -y autoclean && \
    apt-get -y autoremove

RUN mkdir /samtools
WORKDIR /samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2 && \
    bzip2 -d samtools-1.19.2.tar.bz2 && \
    tar -xvf samtools-1.19.2.tar && \
    cd samtools-1.19.2 && \
    autoheader && \
    autoconf -Wno-syntax && \
    ./configure && \
    make && \
    make install

ADD . /kage-lite
WORKDIR /kage-lite

RUN pip install numpy==1.23.5 scipy==1.10.0 dill==0.3.6 pandas==2.1.1 truvari==4.1.0 scikit-learn==1.3.1 matplotlib==3.8.0 pytest==8.0.0 \
                sharedarray==3.2.4 bionumpy==1.0.4 biopython==1.83 npstructures==0.2.16 cython==3.0.8 scikit-allel==1.3.7
RUN pip install -e shared_memory_wrapper/ -e obgraph/ -e graph_kmer_index/ -e kmer_mapper/ -e kage/ -e lite_utils/

ENV NPY_DISABLE_CPU_FEATURES="AVX512F,AVX512_KNL,AVX512_KNM,AVX512_CLX,AVX512_CNL,AVX512_ICL,AVX512CD,AVX512_SKX"

CMD ["bash"]
