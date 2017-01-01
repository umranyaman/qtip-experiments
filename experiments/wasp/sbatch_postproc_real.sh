#!/bin/sh

if [ "$1" = "scavenger" ] ; then
PART1="#SBATCH --partition=scavenger"
PART2="#SBATCH --qos=scavenger"
else
PART1="#SBATCH --partition=shared"
PART2=""
fi

NTHREADS=8

for dat in ERR050082_1 ERR050083_1 ; do
for pe in unp pair ; do
for ext in bwa bt2 snap ; do

MEMGB=6
[ "${ext}" = "bwa" ] && MEMGB=12
[ "${ext}" = "snap" ] && MEMGB=50

P="${dat}.${ext}.${pe}"
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

if [ -f "${P}.wasp_out/${P}.sorted.remap.fq2.gz" ] ; then
    FASTQ1="--fastq  ${P}.wasp_out/${P}.sorted.remap.fq1.gz"
    FASTQ2="--fastq2 ${P}.wasp_out/${P}.sorted.remap.fq2.gz"
fi

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
