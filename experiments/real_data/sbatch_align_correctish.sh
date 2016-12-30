#!/bin/sh

# Do Bowtie 2, BWA and SNAP unpaired runs again but using the "correctish"
# annotated reads

if [ "$1" = "scavenger" ] ; then
PART1="#SBATCH --partition=scavenger"
PART2="#SBATCH --qos=scavenger"
else
PART1="#SBATCH --partition=shared"
PART2=""
fi

for dat in ERR050082_1 ERR050083_1 ; do
for pe in unp ; do

#
# Bowtie 2 and BWA runs
#

for ext in bwa.${pe}.sam bt2.${pe}.sam ; do

cat >.new_${dat}.${ext} <<EOF
#!/bin/sh
#SBATCH --job-name=al_cor_${dat}_${ext}
#SBATCH --output=al_cor_${dat}_${ext}.out
#SBATCH --error=al_cor_${dat}_${ext}.err
${PART1}
${PART2}
#SBATCH --time=8:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=16
/usr/bin/time -v make new_${dat}.${ext}
EOF
echo "sbatch .new_${dat}.${ext}"

done

#
# SNAP runs
#

for ext in snap.${pe}.sam ; do

cat >.new_${dat}.${ext} <<EOF
#!/bin/sh
#SBATCH --job-name=al_cor_${dat}_${ext}
#SBATCH --output=al_cor_${dat}_${ext}.out
#SBATCH --error=al_cor_${dat}_${ext}.err
${PART1}
${PART2}
#SBATCH --time=8:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=16
/usr/bin/time -v make new_${dat}.${ext}
EOF
echo "sbatch .new_${dat}.${ext}"

done

done
done

