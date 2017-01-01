#!/bin/bash -l

#SBATCH
#SBATCH --job-name=MakeSnpDir
#SBATCH --output=.MakeSnpDir.out
#SBATCH --error=.MakeSnpDir.err
#SBATCH --time=24:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1

ENSEMBL_REL=85

URL=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_REL}/variation/vcf/homo_sapiens/Homo_sapiens.vcf.gz
DR=human_ensembl85_snps
if [ ! -d "${DR}" ] ; then
    mkdir -p ${DR}
    rm -f snps
    ln -s -f ${DR} snps
    curl ${URL} | gzip -dc | awk '$1 !~ /^#/ && length($4) == 1 && length($5) == 1 && $4 != $5 && $4 != "." && $5 != "." {print $2,$4,$5 > "snps/"$1".snps.txt"}'
    gzip ${DR}/*.txt
fi

URL=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_REL}/variation/vcf/mus_musculus/Mus_musculus.vcf.gz
DR=mouse_ensembl85_snps
if [ ! -d "${DR}" ] ; then
    mkdir -p ${DR}
    rm -f snps
    ln -s -f ${DR} snps
    curl ${URL} | gzip -dc | awk '$1 !~ /^#/ && length($4) == 1 && length($5) == 1 && $4 != $5 && $4 != "." && $5 != "." {print $2,$4,$5 > "snps/"$1".snps.txt"}'
    gzip ${DR}/*.txt
fi

# Ensembl does not have variation data for Zea Mays
