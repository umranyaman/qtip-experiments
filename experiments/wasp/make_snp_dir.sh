#!/bin/bash -l

#SBATCH
#SBATCH --job-name=MakeSnpDir
#SBATCH --time=24:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1

ENSEMBL_REL=85
URL=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_REL}/variation/vcf/homo_sapiens/Homo_sapiens.vcf.gz

mkdir -p snps
curl ${URL} | gzip -dc | awk '$1 !~ /^#/ && length($4) == 1 && length($5) == 1 && $4 != $5 && $4 != "." && $5 != "." {print $2,$4,$5 > "snps/"$1".snps.txt"}'
gzip snps/*.txt
