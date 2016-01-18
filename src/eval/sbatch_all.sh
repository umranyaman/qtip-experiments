#!/bin/sh

for dat in ERR050082_1 ERR050083_1 ; do
for pe in unp pair ; do

for ext in bt2vs.${pe}.sam bt2vsl.${pe}.sam ; do

cat >.${dat}.${ext} <<EOF
#!/bin/sh
#SBATCH --partition=shared
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
#SBATCH --partition=shared
#SBATCH --time=24:00:00
#SBATCH --mem=12G
/usr/bin/time -v make ${dat}.${ext}
EOF
echo "sbatch .${dat}.${ext}"

done

for ext in snap.${pe}.sam ; do

cat >.${dat}.${ext} <<EOF
#!/bin/sh
#SBATCH --partition=shared
#SBATCH --time=12:00:00
#SBATCH --mem=48G
/usr/bin/time -v make ${dat}.${ext}
EOF
echo "sbatch .${dat}.${ext}"

done

cat >.${dat}.${pe} <<EOF
#!/bin/sh
#SBATCH --partition=shared
#SBATCH --time=12:00:00
#SBATCH --mem=4G
/usr/bin/time -v make ${dat}.unp.csv
EOF
echo "sbatch .${dat}.${pe}"
done
done

