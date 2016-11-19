#!/bin/bash -l

#SBATCH
#SBATCH --job-name=MakeSnpDir
#SBATCH --time=24:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1

DBSNP_RELEASE=144
SNP_FILE=snp${DBSNP_RELEASE}Common.txt
UCSC_COMMON_SNP=http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/${SNP_FILE}

mkdir -p snps

cat >.snpize.awk <<EOF
\$11 == "genomic" && \$12 == "single" {
    rf=\$8;
    alt=substr(\$10,1,1);
    if(rf == alt) {
        alt=substr(\$10,3,1);
    }
    print \$3,rf,alt > "snps/"\$2".snps.txt"
}
EOF

curl $UCSC_COMMON_SNP | \
    gzip -dc |
    awk -v FS='\t' -f .snpize.awk

gzip snps/*.txt
