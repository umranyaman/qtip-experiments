#!/bin/sh

NM=ERR194147
FB=$QTIP_EXPERIMENTS_HOME/software/freebayes/freebayes
FB_BASE="${FB} -i -X -f $QTIP_EXPERIMENTS_HOME/experiments/refs/hg38.fa"
VCFLIB_HOME=$QTIP_EXPERIMENTS_HOME/software/vcflib
VCFISECT=$VCFLIB_HOME/vcflib-git/bin/vcfintersect
MAXDEPTH_FACTOR=4
MAXDEPTH=`python -c "from math import *; print(int(round(50 + ${MAXDEPTH_FACTOR} * sqrt(50))))"`

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
        cat >.CallFB.${NM}_${CHR}_${MINMAPQ}.sh <<EOF
#!/bin/bash -l
#SBATCH
#SBATCH --job-name=CallFB
#SBATCH --output=.CallFB.${NM}_${CHR}_${MINMAPQ}.out
#SBATCH --error=.CallFB.${NM}_${CHR}_${MINMAPQ}.err
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --partition=shared
#SBATCH --time=4:00:00
FN="${NM}_input_${CHR}_${MINMAPQ}"
if [ ! -f ERR194147.sam/\${FN}.raw.vcf ] ; then
    ${FB_BASE} -b ERR194147.sam/input.sorted.bam -r ${CHR} ${MINMAPQ_ARG} -v ERR194147.sam/\${FN}.raw.vcf
fi
if [ ! -f ERR194147.sam/\${FN}.cr_filt.vcf ] ; then
    ${VCFISECT} -b cr_${CHR}.bed ERR194147.sam/\${FN}.raw.vcf | \
        gawk '/^#/ {print} match(\$0, /DP=([0-9]+);/, a) {if(a[1] <= ${MAXDEPTH}) {print}}' > ERR194147.sam/\${FN}.cr_filt.vcf
fi

FN="${NM}_final_${CHR}_${MINMAPQ}"
if [ ! -f ERR194147.sam/\${FN}.raw.vcf ] ; then
    ${FB_BASE} -b ERR194147.sam/final.sorted.bam -r ${CHR} ${MINMAPQ_ARG} -v ERR194147.sam/\${FN}.raw.vcf
fi
if [ ! -f ERR194147.sam/\${FN}.cr_filt.vcf ] ; then
    ${VCFISECT} -b cr_${CHR}.bed ERR194147.sam/\${FN}.raw.vcf | \
        gawk '/^#/ {print} match(\$0, /DP=([0-9]+);/, a) {if(a[1] <= ${MAXDEPTH}) {print}}' > ERR194147.sam/\${FN}.cr_filt.vcf
fi
EOF
        echo "sbatch .CallFB.${NM}_${CHR}_${MINMAPQ}.sh"
        [ "$1" = "wet" ] && sbatch .CallFB.${NM}_${CHR}_${MINMAPQ}.sh && sleep 1
	done
done
