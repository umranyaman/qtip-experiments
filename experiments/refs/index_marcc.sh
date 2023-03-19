#!/bin/sh

[ -z "$QTIP_HOME" ] && echo "Set QTIP_HOME" && exit 1
[ -z "$QTIP_EXPERIMENTS_HOME" ] && echo "Set QTIP_EXPERIMENTS_HOME" && exit 1

do_genome() {
    cat > ".${1}.fa2b.sh" <<EOF
#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -l h_rt=2:00:00
#$ -o .${1}.fa2b.sh.o
#$ -e .${1}.fa2b.sh.e
samtools faidx ${1}.fa
faToTwoBit ${1}.fa ${1}.2bit
EOF
    echo qsub .${1}.fa2b.sh

    cat > .${1}.bt2.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=16G
#$ -l h_rt=10:00:00
#$ -o .${1}.bt2.sh.o
#$ -e .${1}.bt2.sh.e
$QTIP_HOME/software/bowtie2/bowtie2-build --bmax 537647674 --dcv 1024 $QTIP_EXPERIMENTS_HOME/experiments/refs/${1}.fa ${1}.fa
EOF
    echo qsub .${1}.bt2.sh

    cat > .${1}.bwa.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -l h_rt=10:00:00
#$ -o .${1}.bwa.sh.o
#$ -e .${1}.bwa.sh.e
$QTIP_HOME/software/bwa/bwa index $QTIP_EXPERIMENTS_HOME/experiments/refs/${1}.fa
mv $TS_REFS/${1}.fa.* .
EOF
    echo qsub .${1}.bwa.sh

    cat > .${1}.snap.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=60G
#$ -l h_rt=2:00:00
#$ -o .${1}.snap.sh.o
#$ -e .${1}.snap.sh.e
$QTIP_HOME/software/snap/snap/snap-aligner index $QTIP_EXPERIMENTS_HOME/experiments/refs/${1}.fa ${1}.fa.snap -bSpace
EOF
    echo qsub .${1}.snap.sh
}

do_genome hg19
do_genome hg38
do_genome mm10
do_genome zm_AGPv4
