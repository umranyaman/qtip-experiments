#!/bin/sh

ALIGNER_CPUS=$1
[ -z "${ALIGNER_CPUS}" ] && ALIGNER_CPUS=56

mkdir -p temp
${QTIP_HOME}/src/qtip \
    --ref ${QTIP_EXPERIMENTS_HOME}/experiments/refs/hg38.fa \
    --m1 ERR194147_1.fastq \
    --m2 ERR194147_2.fastq \
    --index ${QTIP_EXPERIMENTS_HOME}/experiments/refs/hg38.fa \
    --bt2-exe ${QTIP_HOME}/software/bowtie2/bowtie2 \
    --keep-intermediates \
    --output-directory ERR194147.sam \
    --write-orig-mapq \
    --write-precise-mapq \
    --temp-directory temp \
    -- -I 0 -X 550 -t -p${ALIGNER_CPUS} --reorder

# awk -v FS='\t' '$9 > 0 && $1 !~ /^@/ && $9 < 10000 {h[int($9/10)] += 1} END {for(d in h) {print d*10,h[d]}}' ERR194147_1.fraglen_pairs.sam | sort -r -n -k1,1

