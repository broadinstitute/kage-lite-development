#!/bin/bash

CROMWELL_JAR=$1
OUTPUT_DIR=$2
DOCKER=$3

mkdir -p $2

RESOURCES_DIR=../resources

# run panel WDL
java -Dconfig.file=$RESOURCES_DIR/local.conf -jar $CROMWELL_JAR run KAGEPanelWithPreprocessing.wdl -i <(sed -e "s|__DOCKER__|$DOCKER|g" KAGEPanelWithPreprocessing-chr1-1Mbp-chr2-1Mbp.json) -m $OUTPUT_DIR/metadata.json

INDEX=$(jq -r '.outputs."KAGEPanelWithPreprocessing.index"' $OUTPUT_DIR/metadata.json)
KMER_INDEX=$(jq -r '.outputs."KAGEPanelWithPreprocessing.kmer_index_only_variants_with_revcomp"' $OUTPUT_DIR/metadata.json)
cp $INDEX $OUTPUT_DIR
cp $KMER_INDEX $OUTPUT_DIR

# run kmer mapping and genotyping (w/o helper model)
docker run --shm-size 4G -v $(readlink -m $OUTPUT_DIR):/kage-lite/test $DOCKER \
  kmer_mapper map -d True -c 100000000 \
                  -i /kage-lite/test/$(basename $KMER_INDEX) \
                  -f /kage-lite/wdl/resources/HG00731.final.chr1-1Mbp-chr2-1Mbp.noN.fasta \
                  -o /kage-lite/test/HG00731.final.chr1-1Mbp-chr2-1Mbp.noN.kmer_counts.npy

docker run --shm-size 4G -v $(readlink -m $OUTPUT_DIR):/kage-lite/test $DOCKER \
  kage genotype -s HG00731 \
                -I true \
                --average-coverage 30 \
                -i /kage-lite/test/$(basename $INDEX) \
                -c /kage-lite/test/HG00731.final.chr1-1Mbp-chr2-1Mbp.noN.kmer_counts.npy \
                -o /kage-lite/test/HG00731.final.chr1-1Mbp-chr2-1Mbp.noN.vcf

diff -s $RESOURCES_DIR/HG00731.final.chr1-1Mbp-chr2-1Mbp.noN.expected.vcf $OUTPUT_DIR/HG00731.final.chr1-1Mbp-chr2-1Mbp.noN.vcf