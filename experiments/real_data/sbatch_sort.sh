#!/bin/sh

# run after sbatch_multialign.sh

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

cat >.${dat}.${ext}.${pe}.sort.sh <<EOF
#!/bin/sh
${PART1}
${PART2}
#SBATCH --time=8:00:00
#SBATCH --mem=4G
if samtools sort -o ${dat}.${ext}.${pe}.sorted.bam ${dat}.${ext}.${pe}.sam.bam ; then
    rm -f ${dat}.${ext}.${pe}.sam.bam
fi
EOF
echo "sbatch .${dat}.${ext}.${pe}.sort.sh"

done
done
done
