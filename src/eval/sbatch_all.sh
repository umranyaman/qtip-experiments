#!/bin/sh

for dat in ERR050082_1 ERR050083_1 ; do
for pe in unp pair ; do
for ext in bt2vs.${pe}.sam bt2vsl.${pe}.sam bwa.${pe}.sam bt2.${pe}.sam ; do

cat >.${dat}.${ext} <<EOF
#!/bin/sh
#SBATCH --partition=shared
#SBATCH --time=48:00:00
#SBATCH --mem=12G
/usr/bin/time -v make ${dat}.${ext}
EOF
echo "sbatch .${dat}.${ext}"

: <<'END'
cat >.${dat}.${ext}.plain <<EOF
#!/bin/sh
#SBATCH --partition=shared
#SBATCH --time=48:00:00
#SBATCH --mem=12G
/usr/bin/time -v make ${dat}.${ext}.plain
EOF
echo "sbatch .${dat}.${ext}.plain"
END

done
for ext in snap.${pe}.sam ; do

cat >.${dat}.${ext} <<EOF
#!/bin/sh
#SBATCH --partition=shared
#SBATCH --time=48:00:00
#SBATCH --mem=48G
/usr/bin/time -v make ${dat}.${ext}
EOF
echo "sbatch .${dat}.${ext}"

: <<'END'
cat >.${dat}.${ext}.plain <<EOF
#!/bin/sh
#SBATCH --partition=shared
#SBATCH --time=48:00:00
#SBATCH --mem=48G
/usr/bin/time -v make ${dat}.${ext}.plain
EOF
echo "sbatch .${dat}.${ext}.plain"
END

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
