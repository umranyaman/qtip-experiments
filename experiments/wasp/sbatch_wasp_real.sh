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

SNPS=human_ensembl85_snps

for dat in ERR050082_1 ERR050083_1 ; do
for pe in unp pair ; do
for ext in bwa bt2 snap ; do

PE_ARG=
[ "${pe}" = "pair" ] && PE_ARG="--is_paired_end"

FN="${dat}.${ext}.${pe}"
cat >.${FN}.wasp.sh <<EOF
#!/bin/sh

#SBATCH
#SBATCH --job-name=Wasp
#SBATCH --output Wasp.${FN}.out
#SBATCH --error Wasp.${FN}.err
#SBATCH --nodes=1
${PART1}
${PART2}
#SBATCH --time=8:00:00
#SBATCH --mem=4G

mkdir -p ${FN}.wasp_out
${FIND_SNPS} --is_sorted ${PE_ARG} --output_dir ${FN}.wasp_out --snp_dir ${SNPS} ${FN}.sorted.bam
EOF
echo "sbatch .${FN}.wasp.sh"

done
done
done
