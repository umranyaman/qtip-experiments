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

for rdlen in 100 250 ; do
for genome in hg mm ; do
for aln in bt2s snap bwamem ; do

    # Unpaired
    FN="r0_${aln}_${genome}_mason_ill_${genome}_${rdlen}.trial0"
    PE_ARG=
    cat >.${FN}.wasp.sh <<EOF
#!/bin/sh

#SBATCH
#SBATCH --job-name=Wasp
#SBATCH --output .Wasp.${FN}.out
#SBATCH --error .Wasp.${FN}.err
#SBATCH --nodes=1
${PART1}
${PART2}
#SBATCH --time=8:00:00
#SBATCH --mem=4G

mkdir -p ${FN}.wasp_out
${FIND_SNPS} --is_sorted ${PE_ARG} --output_dir ${FN}.wasp_out --snp_dir ${SNPS} ${FN}.sorted.bam
EOF
    echo "sbatch .${FN}.wasp.sh"

    # Paired
    FN="r12_${aln}${rdlen}_${genome}_mason_ill_${genome}_${rdlen}.trial0"
    PE_ARG="--is_paired_end"
    cat >.${FN}.wasp.sh <<EOF
#!/bin/sh

#SBATCH
#SBATCH --job-name=Wasp
#SBATCH --output .Wasp.${FN}.out
#SBATCH --error .Wasp.${FN}.err
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
