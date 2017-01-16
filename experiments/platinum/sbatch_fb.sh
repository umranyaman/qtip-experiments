#!/bin/sh

NM=ERR194147
FB=$QTIP_EXPERIMENTS_HOME/software/freebayes/freebayes
FB_BASE="${FB} -X -u --haplotype-length 0 -f $QTIP_EXPERIMENTS_HOME/experiments/refs/hg38.fa"
VCFLIB_HOME=$QTIP_EXPERIMENTS_HOME/software/vcflib
VCFISECT=$VCFLIB_HOME/vcflib-git/bin/vcfintersect

for COV in 50 40 30 ; do

    MAXDEPTH_FACTOR=4
    MAXDEPTH=`python -c "from math import *; print(int(round(${COV} + ${MAXDEPTH_FACTOR} * sqrt(${COV}))))"`
    SUBSAMP_FRAC=`python -c "print(${COV}/52.78920049911709)"`

    for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ; do
        for MINMAPQ in 00 01 02 03 04 05 06 07 08 09 10 11 12 15 20 30 d s u ; do
            MINMAPQ_ARG="--min-mapping-quality ${MINMAPQ}"
            if [ "${MINMAPQ}" = "d" ] ; then
                MINMAPQ_ARG=""
            elif [ "${MINMAPQ}" = "s" ] ; then
                MINMAPQ_ARG="--standard-filters"
            elif [ "${MINMAPQ}" = "u" ] ; then
                MINMAPQ_ARG="--use-mapping-quality"
            fi
            LAB="${NM}_${CHR}_${MINMAPQ}_${COV}"
            INP_FN="${NM}.sam/${NM}_input_${CHR}_${MINMAPQ}_${COV}"
            FIN_FN="${NM}.sam/${NM}_final_${CHR}_${MINMAPQ}_${COV}"
            if [ ! -f "${INP_FN}.cr_filt.vcf" -o ! -f "${FIN_FN}.cr_filt.vcf" ] ; then
                cat >.CallFB.${LAB}.sh <<EOF
#!/bin/bash -l
#SBATCH
#SBATCH --job-name=CallFB
#SBATCH --output=.CallFB.${LAB}.out
#SBATCH --error=.CallFB.${LAB}.err
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --partition=shared
#SBATCH --time=8:00:00
if [ ! -f ${INP_FN}.raw.vcf ] ; then
    samtools view -s ${SUBSAMP_FRAC} -bu ERR194147.sam/input.sorted.bam ${CHR} | ${FB_BASE} --stdin ${MINMAPQ_ARG} -v ${INP_FN}.raw.vcf
fi
if [ ! -f ${INP_FN}.cr_filt.vcf ] ; then
    ${VCFISECT} -b cr_${CHR}.bed ${INP_FN}.raw.vcf | \
        gawk '/^#/ {print} match(\$0, /DP=([0-9]+);/, a) {if(a[1] <= ${MAXDEPTH}) {print}}' > ${INP_FN}.cr_filt.vcf
fi

if [ ! -f ${FIN_FN}.raw.vcf ] ; then
    samtools view -s ${SUBSAMP_FRAC} -bu ERR194147.sam/final.sorted.bam ${CHR} | ${FB_BASE} --stdin ${MINMAPQ_ARG} -v ${FIN_FN}.raw.vcf
fi
if [ ! -f ${FIN_FN}.cr_filt.vcf ] ; then
    ${VCFISECT} -b cr_${CHR}.bed ${FIN_FN}.raw.vcf | \
        gawk '/^#/ {print} match(\$0, /DP=([0-9]+);/, a) {if(a[1] <= ${MAXDEPTH}) {print}}' > ${FIN_FN}.cr_filt.vcf
fi
EOF
                echo "sbatch .CallFB.${LAB}.sh"
                [ "$1" = "wet" ] && sbatch .CallFB.${LAB}.sh && sleep 1
            fi
        done
    done
done
