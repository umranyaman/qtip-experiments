#!/bin/sh

BT2_HOME=$HOME/git/bowtie2/bowtie2

REF=$HOME/fasta/hg19/chr1.fa
IDX=$HOME/bowtie2_indexes/Homo_sapiens.GRCh37.60.dna.fa
UNP_IN=$HOME/reads/SRA_HISEQ2000_FC1.shuffle.2M.1.fastq
SAM_OUT=test_ts_human.sam

python ../bin/ts.py \
	--ref $REF \
	--U $UNP_IN \
	--S $SAM_OUT \
	--fastq \
	--scoring "1,2,6,1,5,3,5,3" \
	--num-reads 100 \
	--bt2-exe $BT2_HOME/bowtie2 \
	--bt2-args "-x $IDX" \
	--upto 10000
