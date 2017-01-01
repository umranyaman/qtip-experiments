#!/bin/sh

if [ "$1" = "scavenger" ] ; then
PART1="#SBATCH --partition=scavenger"
PART2="#SBATCH --qos=scavenger"
else
PART1="#SBATCH --partition=shared"
PART2=""
fi

NTHREADS=8

for rdlen in 100 250 ; do
for genome in hg mm ; do
for aln in bt2s snap bwamem ; do

    MEMGB=6
    [ "${aln}" = "bwa" ] && MEMGB=12
    [ "${aln}" = "snap" ] && MEMGB=50

    # Unpaired
    P="r0_${aln}_${genome}_mason_ill_${genome}_${rdlen}.trial0"
    cat >.${P}.postproc.sh <<EOF
#!/bin/sh

#SBATCH
#SBATCH --job-name=Postprocess
#SBATCH --output=.Postprocess.${P}.out
#SBATCH --error=.Postprocess.${P}.err
#SBATCH --nodes=1
${PART1}
${PART2}
#SBATCH --time=4:00:00
#SBATCH --mem=${MEMGB}G
#SBATCH --cpus-per-task=${NTHREADS}

FASTQ1="--fastq ${P}.wasp_out/${P}.sorted.remap.fq.gz"
FASTQ2=""

python postprocess.py \
    --bam ${P}.sorted.bam \
    \${FASTQ1} \${FASTQ2} \
    --threads ${NTHREADS} \
    --output ${P}.csv
EOF
    echo "sbatch .${P}.postproc.sh"

    P="r12_${aln}${rdlen}_${genome}_mason_ill_${genome}_${rdlen}.trial0"
    cat >.${P}.postproc.sh <<EOF
#!/bin/sh

#SBATCH
#SBATCH --job-name=Postprocess
#SBATCH --output=.Postprocess.${P}.out
#SBATCH --error=.Postprocess.${P}.err
#SBATCH --nodes=1
${PART1}
${PART2}
#SBATCH --time=4:00:00
#SBATCH --mem=${MEMGB}G
#SBATCH --cpus-per-task=${NTHREADS}

FASTQ1="--fastq  ${P}.wasp_out/${P}.sorted.remap.fq1.gz"
FASTQ2="--fastq2 ${P}.wasp_out/${P}.sorted.remap.fq2.gz"

python postprocess.py \
    --bam ${P}.sorted.bam \
    \${FASTQ1} \${FASTQ2} \
    --threads ${NTHREADS} \
    --output ${P}.csv
EOF
    echo "sbatch .${P}.postproc.sh"

done
done
done
