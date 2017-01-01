#!/usr/bin/env bash

NTHREADS=8

if [ "$1" = "scavenger" ] ; then
PART1="#SBATCH --partition=scavenger"
PART2="#SBATCH --qos=scavenger"
else
PART1="#SBATCH --partition=shared"
PART2=""
fi

for pa in unp pair ; do
    for aln in bt2 bwa snap ; do
        for dat in 2 3 ; do
            FN="ERR05008${dat}_1.${aln}.${pa}"
            cat >.${FN}.sort.sh <<EOF
#!/bin/sh
#SBATCH --job-name=SortPreWasp
#SBATCH --output SortPreWasp.${FN}.out
#SBATCH --error SortPreWasp.${FN}.err
#SBATCH --cpus-per-task=${NTHREADS}
${PART1}
${PART2}
#SBATCH --time=8:00:00
#SBATCH --mem=4G
if [ ! -f "${FN}.bam" ] ; then
    sambamba view -t ${NTHREADS} -S -f bam ${FN}.sam > ${FN}.bam
fi
if [ ! -f "${FN}.sorted.bam" ] ; then
    if sambamba sort -m 3G -t ${NTHREADS} -o ${FN}.sorted.bam ${FN}.bam ; then
        rm -f ${FN}.bam
    fi
fi
EOF
            echo "sbatch .${FN}.sort.sh"

        done
    done
done

for rdlen in 100 250 ; do
    for genome in hg mm ; do
        for aln in bt2s snap bwamem ; do
            # Unpaired
            FN="r0_${aln}_${genome}_mason_ill_${genome}_${rdlen}.trial0"
            cat >.${FN}.sort.sh <<EOF
#!/bin/sh
#SBATCH --job-name=SortPreWasp
#SBATCH --output SortPreWasp.${FN}.out
#SBATCH --error SortPreWasp.${FN}.err
#SBATCH --cpus-per-task=${NTHREADS}
${PART1}
${PART2}
#SBATCH --time=8:00:00
#SBATCH --mem=4G
if [ ! -f "${FN}.bam" ] ; then
    sambamba view -t ${NTHREADS} -S -f bam ${FN}.sam > ${FN}.bam
fi
if [ ! -f "${FN}.sorted.bam" ] ; then
    if sambamba sort -m 3G -t ${NTHREADS} -o ${FN}.sorted.bam ${FN}.bam ; then
        rm -f ${FN}.bam
    fi
fi
EOF
            echo "sbatch .${FN}.sort.sh"

            # Paired
            FN="r12_${aln}${rdlen}_${genome}_mason_ill_${genome}_${rdlen}.trial0"
            cat >.${FN}.sort.sh <<EOF
#!/bin/sh
#SBATCH --job-name=SortPreWasp
#SBATCH --output SortPreWasp.${FN}.out
#SBATCH --error SortPreWasp.${FN}.err
#SBATCH --cpus-per-task=${NTHREADS}
${PART1}
${PART2}
#SBATCH --time=8:00:00
#SBATCH --mem=4G
if [ ! -f "${FN}.bam" ] ; then
    sambamba view -t ${NTHREADS} -S -f bam ${FN}.sam > ${FN}.bam
fi
if [ ! -f "${FN}.sorted.bam" ] ; then
    if sambamba sort -m 3G -t ${NTHREADS} -o ${FN}.sorted.bam ${FN}.bam ; then
        rm -f ${FN}.bam
    fi
fi
EOF
            echo "sbatch .${FN}.sort.sh"
        done
    done
done
