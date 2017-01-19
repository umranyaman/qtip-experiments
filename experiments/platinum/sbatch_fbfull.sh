#!/bin/sh

NM=ERR194147
FB=$QTIP_EXPERIMENTS_HOME/software/freebayes/freebayes
FB_BASE="${FB} -X -u --haplotype-length 0 -f $QTIP_EXPERIMENTS_HOME/experiments/refs/hg38.fa"
VCFLIB_HOME=$QTIP_EXPERIMENTS_HOME/software/vcflib
VCFISECT=$VCFLIB_HOME/vcflib-git/bin/vcfintersect

#for COV in F 50 40 30 ; do
for COV in F ; do

    MAXDEPTH_FACTOR=4
    SUBSAMP_FRAC=""
    COVNUM=${COV}
    if [ "${COV}" != "F" ] ; then
        SUBSAMP_FRAC=`python -c "print(${COV}/52.78920049911709)"`
        SUBSAMP_FRAC="-s ${SUBSAMP_FRAC}"
    else
        COVNUM=52.78920049911709
    fi
    MAXDEPTH=`python -c "from math import *; print(int(round(${COVNUM} + ${MAXDEPTH_FACTOR} * sqrt(${COVNUM}))))"`

    for MINMAPQ in 00 01 02 03 04 05 06 07 08 09 10 11 12 15 20 30 d s u ; do
        MINMAPQ_ARG="--min-mapping-quality ${MINMAPQ}"
        if [ "${MINMAPQ}" = "d" ] ; then
            MINMAPQ_ARG=""
        elif [ "${MINMAPQ}" = "s" ] ; then
            MINMAPQ_ARG="--standard-filters"
        elif [ "${MINMAPQ}" = "u" ] ; then
            MINMAPQ_ARG="--use-mapping-quality"
        fi
        LAB="${NM}_${MINMAPQ}_${COV}"
        INP_FN="${NM}.sam/${NM}_input_W_${MINMAPQ}_${COV}"
        FIN_FN="${NM}.sam/${NM}_final_W_${MINMAPQ}_${COV}"
        if [ ! -f "${INP_FN}.cr_filt.vcf" -o ! -f "${FIN_FN}.cr_filt.vcf" ] ; then
            cat >.CallFBWhole.${LAB}.sh <<EOF
#!/bin/bash -l
#SBATCH
#SBATCH --job-name=CallFBWhole
#SBATCH --output=.CallFBWhole.${LAB}.out
#SBATCH --error=.CallFBWhole.${LAB}.err
#SBATCH --nodes=1
#SBATCH --mem=32G
#SBATCH --partition=shared
#SBATCH --time=80:00:00
if [ ! -f ${INP_FN}.raw.vcf ] ; then
    ${FB_BASE} ${MINMAPQ_ARG} -v ${INP_FN}.raw.vcf ERR194147.sam/input.sorted.bam
fi
if [ ! -f ${INP_FN}.cr_filt.vcf ] ; then
    ${VCFISECT} -b cr_W.bed ${INP_FN}.raw.vcf | \
        gawk '/^#/ {print} match(\$0, /DP=([0-9]+);/, a) {if(a[1] <= ${MAXDEPTH}) {print}}' > ${INP_FN}.cr_filt.vcf
fi

if [ ! -f ${FIN_FN}.raw.vcf ] ; then
    ${FB_BASE} ${MINMAPQ_ARG} -v ${FIN_FN}.raw.vcf ERR194147.sam/final.sorted.bam
fi
if [ ! -f ${FIN_FN}.cr_filt.vcf ] ; then
    ${VCFISECT} -b cr_W.bed ${FIN_FN}.raw.vcf | \
        gawk '/^#/ {print} match(\$0, /DP=([0-9]+);/, a) {if(a[1] <= ${MAXDEPTH}) {print}}' > ${FIN_FN}.cr_filt.vcf
fi
EOF
            echo "sbatch .CallFBWhole.${LAB}.sh"
            [ "$1" = "wet" ] && sbatch .CallFBWhole.${LAB}.sh && sleep 1
        fi
    done
done
