#!/bin/sh

# make_rec_profiles.sh

# Make some recovery profiles by running sens.py on some genomes with different
# read length parameters

SENS_SCR=../../bin/sens.py
SENS="python $SENS_SCR"

sens() {
	FA=$1
	INDEX=$2
	MIN_LEN=$3
	MAX_LEN=$4
	IDFRAC=$5
	NM=$6
	$SENS \
		--verbose \
		--fasta $1 \
		--scoring "1,2,6,1,5,3,5,3" \
		--num-reads 10000 \
		--min-id $IDFRAC \
		--min-len $MIN_LEN \
		--max-len $MAX_LEN \
		--bt2-exe $HOME/git/bowtie2/bowtie2/bowtie2 \
	 	--bt2-args "-x $INDEX" > ${NM}_recovery.tsv 
}


sens $HOME/bowtie2_indexes/Homo_sapiens.GRCh37.60.dna.fa $HOME/bowtie2_indexes/Homo_sapiens.GRCh37.60.dna.fa 50 50 0.9 hs_rd50_id90
sens $HOME/bowtie2_indexes/Homo_sapiens.GRCh37.60.dna.fa $HOME/bowtie2_indexes/Homo_sapiens.GRCh37.60.dna.fa 100 100 0.9 hs_rd100_id90 
sens $HOME/bowtie2_indexes/Homo_sapiens.GRCh37.60.dna.fa $HOME/bowtie2_indexes/Homo_sapiens.GRCh37.60.dna.fa 200 200 0.9 hs_rd200_id90 
