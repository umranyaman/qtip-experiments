#!/bin/sh

GRCh37=ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
GRCm38=ftp://ftp.ensembl.org/pub/release-75/fasta/mus_musculus/dna/Mus_musculus.GRCm38.75.dna.primary_assembly.fa.gz
ZM=ftp://ftp.ensemblgenomes.org/pub/plants/release-21/fasta/zea_mays/dna/Zea_mays.AGPv3.21.dna.genome.fa.gz

wget $GRCh37

wget $GRCm38

wget $ZM
