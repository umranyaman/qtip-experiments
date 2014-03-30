#!/bin/sh

cat > .hg19.sh <<EOF
#PBS -q batch
#PBS -l walltime=10:00:00
#PBS -j n
#PBS -l pmem=8gb
cd $PWD
bowtie2-build $TS_REFS/hg19.fa hg19.fa
EOF
echo qsub .hg19.sh

cat > .mm10.sh <<EOF
#PBS -q batch
#PBS -l walltime=10:00:00
#PBS -j n
#PBS -l pmem=8gb
cd $PWD
bowtie2-build $TS_REFS/mm10.fa mm10.fa
EOF
echo qsub .mm10.sh

cat > .zm_AGPv3.sh <<EOF
#PBS -q batch
#PBS -l walltime=10:00:00
#PBS -j n
#PBS -l pmem=8gb
cd $PWD
bowtie2-build $TS_REFS/zm_AGPv3.fa zm_AGPv3.fa
EOF
echo qsub .zm_AGPv3.sh
