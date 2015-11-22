#!/bin/sh

AMT=50M
for PAIRED in 0 1 ; do
for LEN in 100 250 ; do

cat >.mason_reads_${PAIRED}_${LEN}.sh <<EOF
#!/bin/bash -l
#SBATCH
#SBATCH --nodes=1
#SBATCH --mem=32G
#SBATCH --partition=shared
#SBATCH --time=10:00:00
make r${PAIRED}_mason_ill_${LEN}_${AMT}.fq.gz
EOF

echo "sbatch .mason_reads_${PAIRED}_${LEN}.sh"

done
done
