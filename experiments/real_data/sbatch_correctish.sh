#!/bin/sh

# Does this need sbatch_sort.sh to have run first?

# Run after sbatch_align.sh and before sbatch_align_correctish.sh

if [ "$1" = "scavenger" ] ; then
PART1="#SBATCH --partition=scavenger"
PART2="#SBATCH --qos=scavenger"
else
PART1="#SBATCH --partition=shared"
PART2=""
fi

for dat in ERR050082 ERR050083 ; do

cat >.${dat}.correctish.sh <<EOF
#!/bin/sh
#SBATCH --job-name=correctish
#SBATCH --output=correctish-%A.out
#SBATCH --error=correctish-%A.err
${PART1}
${PART2}
#SBATCH --time=8:00:00
#SBATCH --mem=4G

pypy correctish.py \
    --sam ${dat}_1.bwa.pair.sam ${dat}_1.bt2.pair.sam ${dat}_1.snap.pair.sam \
    --fastq1 ${dat}_1.fastq --fastq2 ${dat}_2.fastq \
    --output1 new_${dat}_1.fastq --output2 new_${dat}_2.fastq

EOF
echo "sbatch .${dat}.correctish.sh"

done
