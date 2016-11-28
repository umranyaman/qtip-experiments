#!/bin/sh

WASP_DIR="../../software/wasp/wasp-github"
FIND_SNPS="python ${WASP_DIR}/mapping/find_intersecting_snps.py"

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

P="${dat}.${ext}.${pe}"
cat >.${P}.postproc.sh <<EOF
#!/bin/sh

#SBATCH
#SBATCH --job-name=Postprocess
#SBATCH --nodes=1
${PART1}
${PART2}
#SBATCH --time=4:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=${NTHREADS}

FASTQ1="--fastq ${P}.sorted.bam.out/${P}.sorted.remap.fq.gz"
FASTQ2=""

if [ -f "${P}.sorted.bam.out/${P}.sorted.remap.fq2.gz" ] ; then
    FASTQ1="--fastq  ${P}.sorted.bam.out/${P}.sorted.remap.fq1.gz"
    FASTQ2="--fastq2 ${P}.sorted.bam.out/${P}.sorted.remap.fq2.gz"
fi

python postprocess.py \
    --bam ../real_data/${P}.sorted.bam \
    ${FASTQ1} ${FASTQ2} \
    --threads ${NTHREADS} \
    --output ${P}.csv
EOF
echo "sbatch .${P}.postproc.sh"

done
done
done
