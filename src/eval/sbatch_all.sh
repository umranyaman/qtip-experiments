#!/bin/sh

for dat in ERR050082_1 ERR050083_1 ; do
for pe in unp pair ; do
for ext in bt2vs.${pe}.sam bt2vsl.${pe}.sam bwa.${pe}.sam bt2.${pe}.sam ; do
cat >.${dat}.${ext} <<EOF
#!/bin/sh
#SBATCH --partition=shared
#SBATCH --time=48:00:00
#SBATCH --mem=12G
time -v make ${dat}.${ext}
EOF
echo "sbatch .${dat}.${ext}"
done
for ext in snap.${pe}.sam ; do
cat >.${dat}.${ext} <<EOF
#!/bin/sh
#SBATCH --partition=shared
#SBATCH --time=48:00:00
#SBATCH --mem=64G
time -v make ${dat}.${ext}
EOF
echo "sbatch .${dat}.${ext}"
done
done

cat >.${dat}.${pe} <<EOF
#!/bin/sh
#SBATCH --partition=shared
#SBATCH --time=12:00:00
#SBATCH --mem=4G
make ${dat}.unp.csv
EOF
echo "sbatch .${dat}.${pe}"
done

