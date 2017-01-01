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
#SBATCH --output ${FN}.out
#SBATCH --error ${FN}.err
#SBATCH --cpus-per-task=${NTHREADS}
${PART1}
${PART2}
#SBATCH --time=8:00:00
#SBATCH --mem=4G
if sambamba sort -m 3G -t ${NTHREADS} -o ${FN} ${FN}.sam ; then
    rm -f ${FN}.sam
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
#SBATCH --output ${FN}.out
#SBATCH --error ${FN}.err
#SBATCH --cpus-per-task=${NTHREADS}
${PART1}
${PART2}
#SBATCH --time=8:00:00
#SBATCH --mem=4G
if sambamba sort -m 3G -t ${NTHREADS} -o ${FN} ${FN}.sam ; then
    rm -f ${FN}.sam
fi
EOF
            echo "sbatch .${FN}.sort.sh"

            # Paired
            FN="r12_${aln}${rdlen}_${genome}_mason_ill_${genome}_${rdlen}.trial0"
            cat >.${FN}.sort.sh <<EOF
#!/bin/sh
#SBATCH --job-name=SortPreWasp
#SBATCH --output ${FN}.out
#SBATCH --error ${FN}.err
#SBATCH --cpus-per-task=${NTHREADS}
${PART1}
${PART2}
#SBATCH --time=8:00:00
#SBATCH --mem=4G
if sambamba sort -m 3G -t ${NTHREADS} -o ${FN} ${FN}.sam ; then
    rm -f ${FN}.sam
fi
EOF
            echo "sbatch .${FN}.sort.sh"
        done
    done
done
