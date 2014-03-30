#!/bin/sh

cat > .hg19.sh <<EOF
#PBS -q batch
#PBS -l walltime=10:00:00
#PBS -j n
#PBS -l pmem=8gb
#PBS -l vmem=8gb
#PBS -l pvmem=8gb
#PBS -l mem=8gb
cd $PWD
ulimit -v 8388608
bowtie2-build --bmax 537647674 --dcv 1024 $TS_REFS/hg19.fa hg19.fa
EOF
echo qsub .hg19.sh

cat > .mm10.sh <<EOF
#PBS -q batch
#PBS -l walltime=10:00:00
#PBS -j n
#PBS -l pmem=8gb
#PBS -l vmem=8gb
#PBS -l pvmem=8gb
#PBS -l mem=8gb
cd $PWD
ulimit -v 8388608
bowtie2-build --bmax 537647674 --dcv 1024 $TS_REFS/mm10.fa mm10.fa
EOF
echo qsub .mm10.sh

cat > .zm_AGPv3.sh <<EOF
#PBS -q batch
#PBS -l walltime=10:00:00
#PBS -j n
#PBS -l pmem=8gb
#PBS -l vmem=8gb
#PBS -l pvmem=8gb
#PBS -l mem=8gb
cd $PWD
ulimit -v 8388608
bowtie2-build --bmax 537647674 --dcv 1024 $TS_REFS/zm_AGPv3.fa zm_AGPv3.fa
EOF
echo qsub .zm_AGPv3.sh
