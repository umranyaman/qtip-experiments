#!/bin/sh

# Download reference genomes and remove contigs too small for Mason

GRCh38_FA=Homo_sapiens.GRCh38.dna.primary_assembly.fa
GRCh38_URL=ftp://ftp.ensembl.org/pub/release-85/fasta/homo_sapiens/dna/$GRCh38_FA.gz

GRCh37_FA=Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
GRCh37_URL=ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/$GRCh37_FA.gz

GRCm38_FA=Mus_musculus.GRCm38.75.dna.primary_assembly.fa
GRCm38_URL=ftp://ftp.ensembl.org/pub/release-75/fasta/mus_musculus/dna/$GRCm38_FA.gz

ZM_FA=Zea_mays.AGPv4.dna.toplevel.fa
ZM_URL=ftp://ftp.ensemblgenomes.org/pub/plants/release-32/fasta/zea_mays/dna/$ZM_FA.gz

get() {
    if [ ! -f "${1}" ] ; then
        rm -f "${1}.gz"
        wget "${2}"
        gzip -dc "${1}.gz" | pypy remove_short.py "${1}.short" > "${1}" 2> "${1}.lengths"
        samtools faidx "${1}"
    fi
    ln -s -f "${1}" "${3}.fa"
    ln -s -f "${1}.fai" "${3}.fa.fai"
}

get "${GRCh37_FA}" "${GRCh37_URL}" "hg19"
get "${GRCh38_FA}" "${GRCh38_URL}" "hg38"
get "${GRCm38_FA}" "${GRCm38_URL}" "mm10"
get "${ZM_FA}" "${ZM_URL}" "zm_AGPv4"
