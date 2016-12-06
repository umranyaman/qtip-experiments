#!/bin/sh

# all of the following have unpaired/paired and ERR050082/ERR050083 versions
# hence the "4 x"

# 4 x bowtie 2 --very-sensitive
# 4 x bowtie 2 --very-sensitive-local
# 4 x bowtie 2
# 4 x bwa-mem
# 4 x SNAP
# (20 jobs total)

# run sbatch_multialign.sh after this

if [ "$1" = "scavenger" ] ; then
PART1="#SBATCH --partition=scavenger"
PART2="#SBATCH --qos=scavenger"
else
PART1="#SBATCH --partition=shared"
PART2=""
fi

for dat in ERR050082_1 ERR050083_1 ; do
for pe in unp pair ; do

#
# Very sensitive runs
#

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

#
# Bowtie 2 and BWA runs
#

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

if [ "${pe}" = "unp" ] ; then
cat >.new_${dat}.${ext} <<EOF
#!/bin/sh
${PART1}
${PART2}
#SBATCH --time=8:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=16
/usr/bin/time -v make new_${dat}.${ext}
EOF
echo "sbatch .new_${dat}.${ext}"
fi

cat >.${dat}.ext_${ext} <<EOF
#!/bin/sh
${PART1}
${PART2}
#SBATCH --time=8:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=16
/usr/bin/time -v make ${dat}.ext_${ext}
EOF
echo "sbatch .${dat}.ext_${ext}"

done

#
# SNAP runs
#

for ext in snap.${pe}.sam ; do

cat >.${dat}.${ext} <<EOF
#!/bin/sh
${PART1}
${PART2}
#SBATCH --time=12:00:00
#SBATCH --mem=60G
/usr/bin/time -v make ${dat}.${ext}
EOF
echo "sbatch .${dat}.${ext}"

done

done
done

