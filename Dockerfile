FROM condaforge/mambaforge:latest as conda

ADD . /kage-lite
WORKDIR /kage-lite

CMD ["bash"]