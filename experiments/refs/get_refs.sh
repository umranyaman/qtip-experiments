#!/bin/sh

# Download reference genomes and remove contigs too small for Mason

GRCh38_FA=Homo_sapiens.GRCh38.dna.primary_assembly.fa
GRCh38_URL=ftp://ftp.ensembl.org/pub/release-85/fasta/homo_sapiens/dna/$GRCh38_FA.gz

GRCm38_FA=Mus_musculus.GRCm38.75.dna.primary_assembly.fa
GRCm38_URL=ftp://ftp.ensembl.org/pub/release-75/fasta/mus_musculus/dna/$GRCm38_FA.gz

# Note: there appears to be an AGPv4 assembly, but it is not yet published, has
# a Toronto-agreement-style embargo, and the Ensembl website says "The
# underlying assembly and final set of annotations may change until the data
# has been accepted by GenBank."  So I'm laying off it for now.

# http://plants.ensembl.org/Zea_mays/Info/Index

ZM_FA=Zea_mays.AGPv4.dna.toplevel.fa
ZM_URL=ftp://ftp.ensemblgenomes.org/pub/plants/release-32/fasta/zea_mays/dna/$ZM_FA.gz

if [ ! -f $GRCh38_FA ] ; then
    rm -f $GRCh38_FA.gz
    wget $GRCh38_URL
    gzip -dc $GRCh38_FA.gz | python remove_short.py $GRCh38_FA.short > $GRCh38_FA
    samtools faidx $GRCh38_FA
fi
ln -s -f $GRCh38_FA hg38.fa
ln -s -f $GRCh38_FA.fai hg38.fa.fai

if [ ! -f $GRCm38_FA ] ; then
    rm -f $GRCm38_FA.gz
    wget $GRCm38_URL
    gzip -dc $GRCm38_FA.gz | python remove_short.py $GRCm38_FA.short > $GRCm38_FA
    samtools faidx $GRCm38_FA
fi
ln -s -f $GRCm38_FA mm10.fa
ln -s -f $GRCm38_FA.fai mm10.fa.fai

if [ ! -f $ZM_FA ] ; then
    rm -f $ZM_FA.gz
    wget $ZM_URL
    gzip -dc $ZM_FA.gz | python remove_short.py $ZM_FA.short > $ZM_FA
    samtools faidx $ZM_FA
fi
ln -s -f $ZM_FA zm_AGPv4.fa
ln -s -f $ZM_FA.fai zm_AGPv4.fa.fai
