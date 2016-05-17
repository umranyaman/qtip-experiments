#!/bin/sh

# 4 x multi_aligner.py csv file

# sbatch_align.sh jobs should all be completed before running this

if [ "$1" = "scavenger" ] ; then
PART1="#SBATCH --partition=scavenger"
PART2="#SBATCH --qos=scavenger"
else
PART1="#SBATCH --partition=shared"
PART2=""
fi

for dat in ERR050082_1 ERR050083_1 ; do
for pe in unp pair ; do

cat >.${dat}.${pe} <<EOF
#!/bin/sh
${PART1}
${PART2}
#SBATCH --time=12:00:00
#SBATCH --mem=4G
/usr/bin/time -v make ${dat}.${pe}.csv
EOF

echo "sbatch .${dat}.${pe}"

done
done

