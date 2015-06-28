#!/bin/sh

[ -z "$TS_HOME" ] && echo "Set TS_HOME" && exit 1

#
# hg19
#

cat > .hg19.bt2.sh <<EOF
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
echo qsub .hg19.bt2.sh

cat > .hg19.bwa.sh <<EOF
#PBS -q batch
#PBS -l walltime=10:00:00
#PBS -j n
#PBS -l pmem=8gb
#PBS -l vmem=8gb
#PBS -l pvmem=8gb
#PBS -l mem=8gb
cd $PWD
bwa index $TS_REFS/hg19.fa
mv $TS_REFS/hg19.fa.* .
EOF
echo qsub .hg19.bwa.sh

cat > .hg19.mosaik.sh <<EOF
#PBS -q batch
#PBS -l walltime=10:00:00
#PBS -j n
#PBS -l pmem=8gb
#PBS -l vmem=8gb
#PBS -l pvmem=8gb
#PBS -l mem=8gb
cd $PWD
$TS_HOME/software/mosaik/MosaikBuild -fr $TS_REFS/hg19.fa -oa hg19.fa.mosaik
EOF
echo qsub .hg19.mosaik.sh

cat > .hg19.snap.sh <<EOF
#PBS -q batch
#PBS -l walltime=1:00:00
#PBS -j n
#PBS -l pmem=32gb
#PBS -l vmem=32gb
#PBS -l pvmem=32gb
#PBS -l mem=32gb
cd $PWD
$TS_HOME/software/snap/snap/snap-aligner index $TS_REFS/hg19.fa hg19.fa.snap
EOF
echo qsub .hg19.snap.sh

#
# mm10
#

cat > .mm10.bt2.sh <<EOF
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
echo qsub .mm10.bt2.sh

cat > .mm10.bwa.sh <<EOF
#PBS -q batch
#PBS -l walltime=10:00:00
#PBS -j n
#PBS -l pmem=8gb
#PBS -l vmem=8gb
#PBS -l pvmem=8gb
#PBS -l mem=8gb
cd $PWD
bwa index $TS_REFS/mm10.fa
mv $TS_REFS/mm10.fa.* .
EOF
echo qsub .mm10.bwa.sh

cat > .mm10.mosaik.sh <<EOF
#PBS -q batch
#PBS -l walltime=10:00:00
#PBS -j n
#PBS -l pmem=8gb
#PBS -l vmem=8gb
#PBS -l pvmem=8gb
#PBS -l mem=8gb
cd $PWD
$TS_HOME/software/mosaik/MosaikBuild -fr $TS_REFS/mm10.fa -oa mm10.fa.mosaik
EOF
echo qsub .mm10.mosaik.sh

cat > .mm10.snap.sh <<EOF
#PBS -q batch
#PBS -l walltime=1:00:00
#PBS -j n
#PBS -l pmem=32gb
#PBS -l vmem=32gb
#PBS -l pvmem=32gb
#PBS -l mem=32gb
cd $PWD
$TS_HOME/software/snap/snap/snap-aligner index $TS_REFS/mm10.fa mm10.fa.snap
EOF
echo qsub .mm10.snap.sh

#
# zm_AGPv3
#

cat > .zm_AGPv3.bt2.sh <<EOF
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
echo qsub .zm_AGPv3.bt2.sh

cat > .zm_AGPv3.bwa.sh <<EOF
#PBS -q batch
#PBS -l walltime=10:00:00
#PBS -j n
#PBS -l pmem=8gb
#PBS -l vmem=8gb
#PBS -l pvmem=8gb
#PBS -l mem=8gb
cd $PWD
bwa index $TS_REFS/zm_AGPv3.fa
mv $TS_REFS/zm_AGPv3.fa.* .
EOF
echo qsub .zm_AGPv3.bwa.sh

cat > .zm_AGPv3.mosaik.sh <<EOF
#PBS -q batch
#PBS -l walltime=10:00:00
#PBS -j n
#PBS -l pmem=8gb
#PBS -l vmem=8gb
#PBS -l pvmem=8gb
#PBS -l mem=8gb
cd $PWD
$TS_HOME/software/mosaik/MosaikBuild -fr $TS_REFS/zm_AGPv3.fa -oa zm_AGPv3.fa.mosaik
EOF
echo qsub .zm_AGPv3.mosaik.sh

cat > .zm_AGPv3.snap.sh <<EOF
#PBS -q batch
#PBS -l walltime=1:00:00
#PBS -j n
#PBS -l pmem=32gb
#PBS -l vmem=32gb
#PBS -l pvmem=32gb
#PBS -l mem=32gb
cd $PWD
$TS_HOME/software/snap/snap/snap-aligner index $TS_REFS/zm_AGPv3.fa zm_AGPv3.fa.snap
EOF
echo qsub .zm_AGPv3.snap.sh
