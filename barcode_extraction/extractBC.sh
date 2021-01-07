#!/bin/bash

# usage:
# extractBC.sh <basename of reads>

# requires:
# - PANDAseq (version 2.11)
# optional:
# - pigz (version 2.3.4) to speed up file compression
# - sra-toolkit (version 2.10.8) to download read files from SRA


### test if pigz is available and use gzip if not ###
if [ -z "$(command -v pigz)" ]; then
  pigz() {
   gzip $@
  }
  unpigz() {
   gunzip $@
  }
fi


### define variables ###
BCLENGTH=12

# set sample name and directory #
SAMPLE=`basename ${1}`
DIR=`dirname ${1}`


### download reads from SRA ###
if [ ! -f "${DIR}/${SAMPLE}_1.fastq.gz" ]; then
  fastq-dump --split-files --gzip --outdir ${DIR} ${SAMPLE}
fi


### hardtrim barcode reads ###
zcat ${DIR}/${SAMPLE}_1.fastq.gz | awk -v BCLEN=${BCLENGTH} '{if (NR % 2 == 0) print substr($1, 1, BCLEN); else print}' | pigz > ${DIR}/${SAMPLE}_1_trimmed.fq.gz
zcat ${DIR}/${SAMPLE}_2.fastq.gz | awk -v BCLEN=${BCLENGTH} '{if (NR % 2 == 0) print substr($1, 1, BCLEN); else print}' | pigz > ${DIR}/${SAMPLE}_2_trimmed.fq.gz


### assemble reads ###
if [ ! -d "${DIR}/panda_logs" ]; then
  mkdir ${DIR}/panda_logs
fi

pandaseq -f ${DIR}/${SAMPLE}_1_trimmed.fq.gz -r ${DIR}/${SAMPLE}_2_trimmed.fq.gz -d bfsrkm -G ${DIR}/panda_logs/log_${SAMPLE}.txt.bz2 -W ${DIR}/barcodes_${SAMPLE}.fa.bz2


### filter for barcodes of correct length; combine and count unique barcodes ###
bzcat ${DIR}/barcodes_${SAMPLE}.fa.bz2 | awk -v bclen=${BCLENGTH} '{if (NR % 2 == 0 && length($0) == bclen && $1 !~ /N/) print}' | sort | uniq -c | pigz > ${DIR}/barcodes_${SAMPLE}.count.gz
