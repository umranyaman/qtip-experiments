#!/bin/sh

[ -z "$TS_HOME" ] && echo "Set TS_HOME" && exit 1

#
# hg19
#

cat > .hg19.bt2.sh <<EOF
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --partition=shared
#SBATCH --time=10:00:00
#SBATCH --output=.hg19.bt2.sh.o
#SBATCH --error=.hg19.bt2.sh.e
$TS_HOME/software/bowtie2/bowtie2-build --bmax 537647674 --dcv 1024 $TS_REFS/hg19.fa hg19.fa
EOF
echo qsub .hg19.bt2.sh

cat > .hg19.bwa.sh <<EOF
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --partition=shared
#SBATCH --time=10:00:00
#SBATCH --output=.hg19.bwa.sh.o
#SBATCH --error=.hg19.bwa.sh.e
$TS_HOME/software/bwa/bwa index $TS_REFS/hg19.fa
mv $TS_REFS/hg19.fa.* .
EOF
echo qsub .hg19.bwa.sh

cat > .hg19.snap.sh <<EOF
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --partition=shared
#SBATCH --time=2:00:00
#SBATCH --output=.hg19.snap.sh.o
#SBATCH --error=.hg19.snap.sh.e
$TS_HOME/software/snap/snap/snap-aligner index $TS_REFS/hg19.fa hg19.fa.snap -bSpace
EOF
echo qsub .hg19.snap.sh

#
# mm10
#

cat > .mm10.bt2.sh <<EOF
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --partition=shared
#SBATCH --time=10:00:00
#SBATCH --output=.mm10.bt2.sh.o
#SBATCH --error=.mm10.bt2.sh.e
$TS_HOME/software/bowtie2/bowtie2-build --bmax 537647674 --dcv 1024 $TS_REFS/mm10.fa mm10.fa
EOF
echo qsub .mm10.bt2.sh

cat > .mm10.bwa.sh <<EOF
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --partition=shared
#SBATCH --time=10:00:00
#SBATCH --output=.mm10.bwa.sh.o
#SBATCH --error=.mm10.bwa.sh.e
$TS_HOME/software/bwa/bwa index $TS_REFS/mm10.fa
mv $TS_REFS/mm10.fa.* .
EOF
echo qsub .mm10.bwa.sh

cat > .mm10.snap.sh <<EOF
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --partition=shared
#SBATCH --time=2:00:00
#SBATCH --output=.mm10.snap.sh.o
#SBATCH --error=.mm10.snap.sh.e
$TS_HOME/software/snap/snap/snap-aligner index $TS_REFS/mm10.fa mm10.fa.snap -bSpace
EOF
echo qsub .mm10.snap.sh

#
# zm_AGPv3
#

cat > .zm_AGPv3.bt2.sh <<EOF
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --partition=shared
#SBATCH --time=10:00:00
#SBATCH --output=.zm_AGPv3.bt2.sh.o
#SBATCH --error=.zm_AGPv3.bt2.sh.e
$TS_HOME/software/bowtie2/bowtie2-build --bmax 537647674 --dcv 1024 $TS_REFS/zm_AGPv3.fa zm_AGPv3.fa
EOF
echo qsub .zm_AGPv3.bt2.sh

cat > .zm_AGPv3.bwa.sh <<EOF
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --partition=shared
#SBATCH --time=10:00:00
#SBATCH --output=.zm_AGPv3.bwa.sh.o
#SBATCH --error=.zm_AGPv3.bwa.sh.e
$TS_HOME/software/bwa/bwa index $TS_REFS/zm_AGPv3.fa
mv $TS_REFS/zm_AGPv3.fa.* .
EOF
echo qsub .zm_AGPv3.bwa.sh

cat > .zm_AGPv3.snap.sh <<EOF
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --partition=shared
#SBATCH --time=2:00:00
#SBATCH --output=.zm_AGPv3.snap.sh.o
#SBATCH --error=.zm_AGPv3.snap.sh.e
$TS_HOME/software/snap/snap/snap-aligner index $TS_REFS/zm_AGPv3.fa zm_AGPv3.fa.snap -bSpace
EOF
echo qsub .zm_AGPv3.snap.sh
