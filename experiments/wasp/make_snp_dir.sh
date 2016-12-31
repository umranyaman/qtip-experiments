#!/bin/bash -l

#SBATCH
#SBATCH --job-name=MakeSnpDir
#SBATCH --time=24:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1

ENSEMBL_REL=85
URL=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_REL}/variation/vcf/homo_sapiens/Homo_sapiens.vcf.gz

mkdir -p human_ensembl85_snps
curl ${URL} | gzip -dc | awk '$1 !~ /^#/ && length($4) == 1 && length($5) == 1 && $4 != $5 && $4 != "." && $5 != "." {print $2,$4,$5 > "human_ensembl85_snps/"$1".snps.txt"}'
gzip human_ensembl85_snps/*.txt

URL=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_REL}/variation/vcf/mus_musculus/Mus_musculus.vcf.gz

mkdir -p mouse_ensembl85_snps
curl ${URL} | gzip -dc | awk '$1 !~ /^#/ && length($4) == 1 && length($5) == 1 && $4 != $5 && $4 != "." && $5 != "." {print $2,$4,$5 > "mouse_ensembl85_snps/"$1".snps.txt"}'
gzip mouse_ensembl85_snps/*.txt

# Ensembl does not have variation data for Zea Mays
