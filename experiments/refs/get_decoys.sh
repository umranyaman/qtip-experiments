#!/bin/sh

#
# Get GRCh38 decoys
#

URL="https://sourceforge.net/projects/bio-bwa/files/bwakit"
VER="0.7.12"
NM="bwakit-${VER}"
AR="${NM}-linux.tar.bz2"
DL="${URL}/${AR}/download"

if [ ! -f "hg38_decoy.fa" ] ; then
    wget -O "${AR}" "${DL}"
    tar xvfj "${AR}"
    cp bwa.kit/resource-GRCh38/hs38DH-extra.fa hg38_decoy.fa
    rm -rf bwa.kit
fi

#
# Get GRCh37.p13 decoys from 1K genomes
#

URL_1KG="ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence"
AR_1KG="hs37d5ss.fa.gz"

if [ ! -f "hg19_decoy.fa" ] ; then
    wget -O "${AR_1KG}" "${URL_1KG}/${AR_1KG}"
    gzip -dc "${AR_1KG}" > hg19_decoy.fa
    rm -f "${AR_1KG}"
fi
