#!/bin/bash

# usage:
# subassembly.sh <basename of reads> <reference sequence file>

# requires:
# - bowtie2 (version 2.4.1)
# - bioawk (https://github.com/lh3/bioawk; commit fd40150)
# - cutadapt (version 2.5)
# - trim_galore (version 0.6.6)
# - PANDAseq (version 2.11)
# - R (version 4.0.0) with libraries:
#   - dplyr (version 1.0.0)
#   - readr (version 1.3.1)
#   - stringr (version 1.4.0)
# optional:
# - pigz (version 2.3.4) to speed up file compression
# - sra-toolkit (version 2.10.8) to download read files from SRA


### define variables ###
SAMPLE=`basename ${1}`
DIR=`dirname ${1}`

REFSEQ="${2}"
INDEX="${DIR}/bowtie2-index/$(basename ${REFSEQ%.*})"

BCLENGTH=12

NSLOTS=${NSLOTS:-1}


### test if pigz is available and use gzip if not ###
if [ -z "$(command -v pigz)" ]; then
  pigz() {
   gzip $@
  }
  unpigz() {
   gunzip $@
  }
fi


### build bowtie2 index ###
if [ ! -f "${INDEX}.1.bt2" ]; then
  INDEXDIR=`dirname ${INDEX}`

  if [ ! -d "${INDEXDIR}" ]; then
    mkdir -p ${INDEXDIR}
  fi

  bowtie2-build ${REFSEQ} ${INDEX}
fi


### download reads from SRA ###
if [ ! -f "${DIR}/${SAMPLE}_1.fastq.gz" ]; then
  fastq-dump --split-files --gzip --outdir ${DIR} ${SAMPLE}
fi


### trim reads (required for the 35S minimal promoter and promoters truncated during synthesis) ###
trim_galore --quality 20 --length 30 --adapter file:adapters.fa --adapter2 CAGCTCCTGGAGACCACTAGGCGCGCCTTG --basename ${SAMPLE} --output_dir ${DIR} --paired ${DIR}/${SAMPLE}_1.fastq.gz ${DIR}/${SAMPLE}_4.fastq.gz


### assemble reads ###
if [ ! -d "${DIR}/panda_logs" ]; then
  mkdir ${DIR}/panda_logs
fi

pandaseq -f ${DIR}/${SAMPLE}_val_1.fq.gz -r ${DIR}/${SAMPLE}_val_2.fq.gz -d bfsrkm -G ${DIR}/panda_logs/log_promoters_${SAMPLE}.txt.bz2 -W ${DIR}/promoters_${SAMPLE}.fa.bz2


### align assembled reads to promoter array ###
bowtie2 -p ${NSLOTS} --xeq --no-unal -f --norc -x ${INDEX} -U ${DIR}/promoters_${SAMPLE}.fa.bz2 2> ${DIR}/alignment_stats.txt | pigz > ${DIR}/${SAMPLE}.sam.gz


### extract info from sam file ###
zcat ${DIR}/${SAMPLE}.sam.gz | bioawk -c sam -t 'BEGIN{print "read", "gene", "start", "cigar", "sequence"} {sub(/:[ACGTN+]*;.*$/, "", $qname); print $qname, $rname, $pos, $cigar, $seq}' > ${DIR}/promoters_${SAMPLE}.tsv


### hardtrim barcode reads ###
zcat ${DIR}/${SAMPLE}_2.fastq.gz | awk -v BCLEN=${BCLENGTH} '{if (NR % 2 == 0) print substr($1, 1, BCLEN); else print}' | pigz > ${DIR}/barcode_${SAMPLE}_1.fq.gz
zcat ${DIR}/${SAMPLE}_3.fastq.gz | awk -v BCLEN=${BCLENGTH} '{if (NR % 2 == 0) print substr($1, 1, BCLEN); else print}' | pigz > ${DIR}/barcode_${SAMPLE}_2.fq.gz


### assemble barcodes ###
pandaseq -f ${DIR}/barcode_${SAMPLE}_1.fq.gz -r ${DIR}/barcode_${SAMPLE}_2.fq.gz -d bfsrkm -G ${DIR}/panda_logs/log_barcodes_${SAMPLE}.txt.bz2 -W ${DIR}/barcodes_${SAMPLE}.fa.bz2


### filter for barcodes of correct length and without uncalled bases ###
bzcat ${DIR}/barcodes_${SAMPLE}.fa.bz2 | awk -v BCLEN=${BCLENGTH} -v OFS='\t' 'BEGIN{print "read", "barcode"} {sub(/:[ACGTN+]*;.*$/, "", $1); NAME=substr($1, 2); getline; if (length($1) == BCLEN && $1 ~ /^[ACGT]+$/) print NAME, $0}' > ${DIR}/barcodes_${SAMPLE}.tsv


### merge barcode and fragment files + combine and count unique fragments (min 5 reads) ###
Rscript join_bc_seq.R ${DIR}/barcodes_${SAMPLE}.tsv ${DIR}/promoters_${SAMPLE}.tsv ${DIR}/subassembly_${SAMPLE}.tsv 5


### convert alignment info to variant name ###
Rscript cigar_to_variant.R ${DIR}/subassembly_${SAMPLE}.tsv ${REFSEQ}


### compress final files ###
pigz ${DIR}/promoters_${SAMPLE}.tsv ${DIR}/barcodes_${SAMPLE}.tsv ${DIR}/subassembly_${SAMPLE}.tsv ${DIR}/subassembly_${SAMPLE}_variant.tsv