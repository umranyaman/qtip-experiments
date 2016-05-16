#!/bin/sh

# 24 jobs in here

# all of the following have unpaired/paired and ERR050082/ERR050083 versions
# hence the "4 x"

# Alignment jobs:
# 4 x bowtie 2 --very-sensitive
# 4 x bowtie 2 --very-sensitive-local
# 4 x bowtie 2
# 4 x bwa-mem
# 4 x SNAP

# Other:
# 4 x multi_aligner.py csv file

# multi_aligner.py jobs should only be run once all alignment jobs are done

if [ "$1" = "scavenger" ] ; then
PART1="#SBATCH --partition=scavenger"
PART2="#SBATCH --qos=scavenger"
else
PART1="#SBATCH --partition=shared"
PART2=""
fi

for dat in ERR050082_1 ERR050083_1 ; do
for pe in unp pair ; do

for ext in bt2vs.${pe}.sam bt2vsl.${pe}.sam ; do

cat >.${dat}.${ext} <<EOF
#!/bin/sh
${PART1}
${PART2}
#SBATCH --time=8:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=16
/usr/bin/time -v make ${dat}.${ext}
EOF
echo "sbatch .${dat}.${ext}"

done

for ext in bwa.${pe}.sam bt2.${pe}.sam ; do

cat >.${dat}.${ext} <<EOF
#!/bin/sh
${PART1}
${PART2}
#SBATCH --time=24:00:00
#SBATCH --mem=12G
/usr/bin/time -v make ${dat}.${ext}
EOF
echo "sbatch .${dat}.${ext}"

done

for ext in snap.${pe}.sam ; do

cat >.${dat}.${ext} <<EOF
#!/bin/sh
${PART1}
${PART2}
#SBATCH --time=12:00:00
#SBATCH --mem=48G
/usr/bin/time -v make ${dat}.${ext}
EOF
echo "sbatch .${dat}.${ext}"

done

cat >.${dat}.${pe} <<EOF
#!/bin/sh
${PART1}
${PART2}
#SBATCH --time=12:00:00
#SBATCH --mem=4G
/usr/bin/time -v make ${dat}.unp.csv
EOF
echo "sbatch .${dat}.${pe}"
done
done

