#!/bin/sh

GRCh37_FA=Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
GRCh37_URL=ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/$GRCh37_FA.gz

GRCm38_FA=Mus_musculus.GRCm38.75.dna.primary_assembly.fa
GRCm38_URL=ftp://ftp.ensembl.org/pub/release-75/fasta/mus_musculus/dna/$GRCm38_FA.gz

ZM_FA=Zea_mays.AGPv3.21.dna.genome.fa
ZM_URL=ftp://ftp.ensemblgenomes.org/pub/plants/release-21/fasta/zea_mays/dna/$ZM_FA.gz

if [ ! -f $GRCh37_FA ] ; then
    rm -f $GRCh37_FA.gz
    wget $GRCh37_URL
    gunzip $GRCh37_FA.gz
    samtools faidx $GRCh37_FA
fi
ln -s -f $GRCh37_FA hg19.fa

if [ ! -f $GRCm38_FA ] ; then
    rm -f $GRCm38_FA.gz
    wget $GRCm38_URL
    gunzip $GRCm38_FA.gz
    samtools faidx $GRCm38_FA
fi
ln -s -f $GRCm38_FA mm10.fa

if [ ! -f $ZM_FA ] ; then
    rm -f $ZM_FA.gz
    wget $ZM_URL
    gunzip $ZM_FA.gz
    samtools faidx $ZM_FA
fi
ln -s -f $ZM_FA zm_AGPv3.fa
