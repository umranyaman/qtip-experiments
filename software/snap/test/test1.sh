#!/bin/sh

BASE_URL=http://www.cs.jhu.edu/~langmea/resources

[ ! -f lambda_virus.fa ] && wget $BASE_URL/lambda_virus.fa
[ ! -f reads_1.fq ] && wget $BASE_URL/reads_1.fq
[ ! -f reads_2.fq ] && wget $BASE_URL/reads_2.fq

# Index
../snap/snap-aligner index lambda_virus.fa lambda_virus.snap

# Simple unpaired
../snap/snap-aligner single lambda_virus.snap reads_1.fq -o -sam unp1.sam

# Unpaired using stdin
cat reads_1.fq | ../snap/snap-aligner single lambda_virus.snap -fastq - -o -sam unp2.sam

# Paired-end, separate files
../snap/snap-aligner paired lambda_virus.snap reads_1.fq reads_2.fq -o -sam pai1.sam

python ../../../bin/fastq_interleave.py reads_1.fq reads_2.fq > reads_paired.fq

# Paired-end, interleaved file
../snap/snap-aligner paired lambda_virus.snap -pairedInterleavedFastq reads_paired.fq -o -sam pai2.sam

# Paired-end, interleaved, stdin file
cat reads_paired.fq | ../snap/snap-aligner paired lambda_virus.snap -pairedInterleavedFastq - -o -sam pai3.sam

# Paired-end twice
../snap/snap-aligner paired lambda_virus.snap reads_1.fq reads_2.fq reads_1.fq reads_2.fq -o -sam pai_d1.sam
