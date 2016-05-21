#!/bin/sh

# Download genomes that are known in the literature to be contaminants

# Diverse and Widespread Contamination Evident in the Unmapped Depths of High Throughput Sequencing Data
# doi:10.1371/journal.pone.0110808

# Assessing the prevalence of mycoplasma contamination in cell culture via a survey of NCBI’s RNA-seq archive
# doi:10.1093/nar/gkv136

BACT="ftp://ftp.ncbi.nih.gov/genomes/genbank/bacteria/"
FUNGI="ftp://ftp.ncbi.nih.gov/genomes/genbank/fungi/"

# Propionibacterium acnes
FA="GCA_001481615.1_ASM148161v1_genomic.fna"
if [ ! -f "${FA}" ] ; then
    wget "${BACT}/Propionibacterium_acnes/latest_assembly_versions/GCA_001481615.1_ASM148161v1/${FA}.gz"
    gzip -dc "${FA}.gz" > "${FA}"
fi

# Malassezia_globosa
FA="GCA_001264815.1_ASM126481v1_genomic.fna"
if [ ! -f "${FA}" ] ; then
    wget "${FUNGI}/Malassezia_globosa/latest_assembly_versions/GCA_001264815.1_ASM126481v1/${FA}.gz"
    gzip -dc "${FA}.gz" > "${FA}"
fi

# Mycoplasma hominis
FA="GCA_001063305.1_ASM106330v1_genomic.fna"
if [ ! -f "${FA}" ] ; then
    wget "${BACT}/Mycoplasma_hominis/latest_assembly_versions/GCA_001063305.1_ASM106330v1/${FA}.gz"
    gzip -dc "${FA}.gz" > "${FA}"
fi

# Mycoplasma hyorhinis
FA="GCA_000496815.1_ASM49681v1_genomic.fna"
if [ ! -f "${FA}" ] ; then
    wget "${BACT}/Mycoplasma_hyorhinis/latest_assembly_versions/GCA_000496815.1_ASM49681v1/${FA}.gz"
    gzip -dc "${FA}.gz" > "${FA}"
fi

# Mycoplasma fermentans
FA="GCA_000209735.1_ASM20973v1_genomic.fna"
if [ ! -f "${FA}" ] ; then
    wget "${BACT}/Mycoplasma_fermentans/latest_assembly_versions/GCA_000209735.1_ASM20973v1/${FA}.gz"
    gzip -dc "${FA}.gz" > "${FA}"
fi

# Acholeplasma laidlawii
FA="GCA_000018785.1_ASM1878v1_genomic.fna"
if [ ! -f "${FA}" ] ; then
    wget "${BACT}/Acholeplasma_laidlawii/all_assembly_versions/GCA_000018785.1_ASM1878v1/${FA}.gz"
    gzip -dc "${FA}.gz" > "${FA}"
fi
