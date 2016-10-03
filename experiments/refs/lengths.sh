#!/bin/bash

GRCh37_FA=Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
GRCh37_URL=ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/$GRCh37_FA.gz

GRCh38_FA=Homo_sapiens.GRCh38.dna.primary_assembly.fa
GRCh38_URL=ftp://ftp.ensembl.org/pub/release-85/fasta/homo_sapiens/dna/$GRCh38_FA.gz

GRCm38_FA=Mus_musculus.GRCm38.75.dna.primary_assembly.fa
GRCm38_URL=ftp://ftp.ensembl.org/pub/release-75/fasta/mus_musculus/dna/$GRCm38_FA.gz

ZM_FA=Zea_mays.AGPv3.21.dna.genome.fa
ZM_URL=ftp://ftp.ensemblgenomes.org/pub/plants/release-21/fasta/zea_mays/dna/$ZM_FA.gz

ZM4_FA=Zea_mays.AGPv4.dna.toplevel.fa
ZM4_URL=ftp://ftp.ensemblgenomes.org/pub/plants/release-32/fasta/zea_mays/dna/$ZM4_FA.gz

#curl $GRCh37_URL | gzip -dc | grep -v '^>' | tee >(tr -cd acgtACGT | wc -c | tee grch37_acgt_length.txt | xargs echo "GRCh37 ACGT len") | \
#                                                   tr -cd acgtnACGTN | wc -c | tee grch37_acgtn_length.txt | xargs echo "GRCh37 ACGTN len"

curl $GRCh38_URL | gzip -dc | grep -v '^>' | tee >(tr -cd acgtACGT | wc -c | tee grch38_acgt_length.txt | xargs echo "GRCh38 ACGT len") | \
                                                   tr -cd acgtnACGTN | wc -c | tee grch38_acgtn_length.txt | xargs echo "GRCh38 ACGTN len"

curl $GRCm38_URL | gzip -dc | grep -v '^>' | tee >(tr -cd acgtACGT | wc -c | tee grcm38_acgt_length.txt | xargs echo "GRCm38 ACGT len") | \
                                                   tr -cd acgtnACGTN | wc -c | tee grcm38_acgtn_length.txt | xargs echo "GRCm38 ACGTN len"

#curl $ZM_URL     | gzip -dc | grep -v '^>' | tee >(tr -cd acgtACGT | wc -c | tee zm_acgt_length.txt | xargs echo "ZM ACGT len") | \
#                                                   tr -cd acgtnACGTN | wc -c | tee zm_acgtn_length.txt | xargs echo "ZM ACGTN len"

curl $ZM4_URL    | gzip -dc | grep -v '^>' | tee >(tr -cd acgtACGT | wc -c | tee zm4_acgt_length.txt | xargs echo "ZM4 ACGT len") | \
                                                   tr -cd acgtnACGTN | wc -c | tee zm4_acgtn_length.txt | xargs echo "ZM4 ACGTN len"
