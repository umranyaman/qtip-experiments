#!/bin/sh

for dat in ERR050082_1 ERR050083_1 ; do
for ext in bt2vs.pair.sam bt2vsl.pair.sam bwa.pair.sam bt2.pair.sam ; do
cat >.${dat}.${ext} <<EOF
#!/bin/sh
#SBATCH --partition=shared
#SBATCH --time=48:00:00
#SBATCH --mem=12G
time -v make ${dat}.${ext}
EOF
echo "sbatch .${dat}.${ext}"
done
for ext in snap.pair.sam ; do
cat >.${dat}.${ext} <<EOF
#!/bin/sh
#SBATCH --partition=shared
#SBATCH --time=48:00:00
#SBATCH --mem=64G
time -v make ${dat}.${ext}
EOF
echo "sbatch .${dat}.${ext}"
done

cat >.${dat}.${ext} <<EOF
#SBATCH --partition=shared
#SBATCH --time=12:00:00
#SBATCH --mem=4G
make ${dat}.unp.csv
EOF
echo "sbatch .${dat}.${ext}"
done

