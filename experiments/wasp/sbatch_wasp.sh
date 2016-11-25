#!/bin/sh

WASP_VER=0.2.1
WASP_DIR="../../software/wasp/WASP-${WASP_VER}"
FIND_SNPS="python ${WASP_DIR}/mapping/find_intersecting_snps.py"

if [ "$1" = "scavenger" ] ; then
PART1="#SBATCH --partition=scavenger"
PART2="#SBATCH --qos=scavenger"
else
PART1="#SBATCH --partition=shared"
PART2=""
fi

for dat in ERR050082_1 ERR050083_1 ; do
for pe in unp pair ; do
for ext in bwa bt2 snap ; do

cat >.${dat}.${ext}.${pe}.wasp.sh <<EOF
#!/bin/sh

#SBATCH
#SBATCH --job-name=Wasp
#SBATCH --nodes=1
${PART1}
${PART2}
#SBATCH --time=8:00:00
#SBATCH --mem=4G

mkdir -p ${dat}.${ext}.${pe}.sorted.bam.out
$FIND_SNPS --is_sorted --output_dir ${dat}.${ext}.${pe}.sorted.bam.out --snp_dir snps ../real_data/${dat}.${ext}.${pe}.sorted.bam
EOF
echo "sbatch .${dat}.${ext}.${pe}.wasp.sh"

done
done
done
