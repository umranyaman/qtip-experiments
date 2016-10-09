#!/bin/sh

[ -z "$QTIP_HOME" ] && echo "Set QTIP_HOME" && exit 1
[ -z "$QTIP_EXPERIMENTS_HOME" ] && echo "Set QTIP_EXPERIMENTS_HOME" && exit 1

#
# hg38
#

cat > .hg38.fa2b.sh <<EOF
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --partition=shared
#SBATCH --time=2:00:00
#SBATCH --output=.hg38.fa2b.sh.o
#SBATCH --error=.hg38.fa2b.sh.e
samtools faidx hg38.fa
faToTwoBit hg38.fa hg38.2bit
EOF
echo sbatch .hg38.fa2b.sh

cat > .hg38.bt2.sh <<EOF
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --partition=shared
#SBATCH --time=10:00:00
#SBATCH --output=.hg38.bt2.sh.o
#SBATCH --error=.hg38.bt2.sh.e
$QTIP_HOME/software/bowtie2/bowtie2-build --bmax 537647674 --dcv 1024 $QTIP_EXPERIMENTS_HOME/experiments/refs/hg38.fa hg38.fa
EOF
echo sbatch .hg38.bt2.sh

cat > .hg38.bwa.sh <<EOF
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --partition=shared
#SBATCH --time=10:00:00
#SBATCH --output=.hg38.bwa.sh.o
#SBATCH --error=.hg38.bwa.sh.e
$QTIP_HOME/software/bwa/bwa index $QTIP_EXPERIMENTS_HOME/experiments/refs/hg38.fa
mv $TS_REFS/hg38.fa.* .
EOF
echo sbatch .hg38.bwa.sh

cat > .hg38.snap.sh <<EOF
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --partition=shared
#SBATCH --time=2:00:00
#SBATCH --output=.hg38.snap.sh.o
#SBATCH --error=.hg38.snap.sh.e
$QTIP_HOME/software/snap/snap/snap-aligner index $QTIP_EXPERIMENTS_HOME/experiments/refs/hg38.fa hg38.fa.snap -bSpace
EOF
echo sbatch .hg38.snap.sh

#
# mm10
#

cat > .mm10.fa2b.sh <<EOF
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --partition=shared
#SBATCH --time=10:00:00
#SBATCH --output=.mm10.fa2b.sh.o
#SBATCH --error=.mm10.fa2b.sh.e
samtools faidx mm10.fa
faToTwoBit mm10.fa mm10.2bit
EOF
echo sbatch .mm10.fa2b.sh

cat > .mm10.bt2.sh <<EOF
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --partition=shared
#SBATCH --time=10:00:00
#SBATCH --output=.mm10.bt2.sh.o
#SBATCH --error=.mm10.bt2.sh.e
$QTIP_HOME/software/bowtie2/bowtie2-build --bmax 537647674 --dcv 1024 $QTIP_EXPERIMENTS_HOME/experiments/refs/mm10.fa mm10.fa
EOF
echo sbatch .mm10.bt2.sh

cat > .mm10.bwa.sh <<EOF
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --partition=shared
#SBATCH --time=10:00:00
#SBATCH --output=.mm10.bwa.sh.o
#SBATCH --error=.mm10.bwa.sh.e
$QTIP_HOME/software/bwa/bwa index $QTIP_EXPERIMENTS_HOME/experiments/refs/mm10.fa
mv $TS_REFS/mm10.fa.* .
EOF
echo sbatch .mm10.bwa.sh

cat > .mm10.snap.sh <<EOF
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --partition=shared
#SBATCH --time=2:00:00
#SBATCH --output=.mm10.snap.sh.o
#SBATCH --error=.mm10.snap.sh.e
$QTIP_HOME/software/snap/snap/snap-aligner index $QTIP_EXPERIMENTS_HOME/experiments/refs/mm10.fa mm10.fa.snap -bSpace
EOF
echo sbatch .mm10.snap.sh

#
# zm_AGPv4
#

cat > .zm_AGPv4.fa2b.sh <<EOF
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --partition=shared
#SBATCH --time=10:00:00
#SBATCH --output=.zm_AGPv4.fa2b.sh.o
#SBATCH --error=.zm_AGPv4.fa2b.sh.e
samtools faidx zm_AGPv4.fa
faToTwoBit zm_AGPv4.fa zm_AGPv4.2bit
EOF
echo sbatch .zm_AGPv4.fa2b.sh

cat > .zm_AGPv4.bt2.sh <<EOF
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --partition=shared
#SBATCH --time=10:00:00
#SBATCH --output=.zm_AGPv4.bt2.sh.o
#SBATCH --error=.zm_AGPv4.bt2.sh.e
$QTIP_HOME/software/bowtie2/bowtie2-build --bmax 537647674 --dcv 1024 $QTIP_EXPERIMENTS_HOME/experiments/refs/zm_AGPv4.fa zm_AGPv4.fa
EOF
echo sbatch .zm_AGPv4.bt2.sh

cat > .zm_AGPv4.bwa.sh <<EOF
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --partition=shared
#SBATCH --time=10:00:00
#SBATCH --output=.zm_AGPv4.bwa.sh.o
#SBATCH --error=.zm_AGPv4.bwa.sh.e
$QTIP_HOME/software/bwa/bwa index $QTIP_EXPERIMENTS_HOME/experiments/refs/zm_AGPv4.fa
mv $QTIP_EXPERIMENTS_HOME/experiments/refs/zm_AGPv4.fa.* .
EOF
echo sbatch .zm_AGPv4.bwa.sh

cat > .zm_AGPv4.snap.sh <<EOF
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --partition=shared
#SBATCH --time=2:00:00
#SBATCH --output=.zm_AGPv4.snap.sh.o
#SBATCH --error=.zm_AGPv4.snap.sh.e
$QTIP_HOME/software/snap/snap/snap-aligner index $QTIP_EXPERIMENTS_HOME/experiments/refs/zm_AGPv4.fa zm_AGPv4.fa.snap -bSpace
EOF
echo sbatch .zm_AGPv4.snap.sh
