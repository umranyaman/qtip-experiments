#!/bin/sh

#
# Make file with CHM1 "decoy" sequences w/r/t GRCh37/hg19
#

# Get assembly
ASM=GCA_000772585.3_ASM77258v3
ASM_AR=${ASM}_genomic.fna
ARGZ=${ASM_AR}.gz
ASM_URL=http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/genomes/genbank/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/${ASM}/${ARGZ}

if [ ! -f ${ASM_AR} ] ; then
    wget ${ASM_URL}
    gunzip ${ARGZ}
fi

# Make genome file
if [ ! -f "${ASM_AR}.genome" ] ; then
    samtools faidx ${ASM_AR}
    awk -v OFS='\t' '{print $1,$2}' "${ASM_AR}.fai" > "${ASM_AR}.genome"
fi

# Get assemblytics files
AR_37=Homo_sapiens_MHAP_assembly.Assemblytics_results.zip
URL_37=http://assemblytics.com/user_data/human/${AR_37}
BED_37=assemblytics_37/JvA9dvkVNwyE4UQ6Q7sj/Homo_sapiens_MHAP_assembly.Assemblytics_structural_variants.bed

ID_38=sIOvqHUKNaCpU8NXdcTS
AR_38=CHM1htert_GRCh38.Assemblytics_results.zip
URL_38=http://assemblytics.com/user_data/${ID_38}/${AR_38}
BED_38=assemblytics_38/${ID_38}/CHM1htert_GRCh38.Assemblytics_structural_variants.bed

# Unpack
unpack() {
    if [ ! -d $1 ] ; then
        wget $2
        unzip $3
        mv user_data assemblytics_37
    fi
}

unpack assemblytics_37 ${URL_37} ${AR_37}
unpack assemblytics_38 ${URL_38} ${AR_38}

getbed() {
    # Make bed with CHM1 insertions using CHM1 coordinates
    BED=$1
    SLOP=$2
    MIN=$3
    if [ ! -f "${4}" ] ; then
        # This starts with the assemblytics bed and does several things:
        # 1. Removes header
        # 2. Removes non-insertion lines
        # 3. Extracts column 10, with the CHM1 coordinates
        # 4. Parses column 10 and converts it to 3-column bed
        # 5. Filters out intervals shorter than 50
        awk '{n += 1; if(n > 1) {if($7 == "Insertion") {print $10}}}' ${BED} | sed 's/[:-]/,/g' | awk -v FS=',' -v OFS='\t' '$3-$2 >= 50 {print $1,$2,$3}' > .tmp.bed
        # Now add slop on both sides
        bedtools slop -i .tmp.bed -g "${ASM_AR}.genome" -b ${SLOP} > "${4}"
    fi

    # Now extract sequences from CHM1
    if [ ! -f "${5}" ] ; then
        bedtools getfasta -fi ${ASM_AR} -bed "${4}" > "${5}"
    fi
}

getbed ${BED_37} 25 50 chm1_grch19.bed hg19_chm1.fa
getbed ${BED_38} 25 50 chm1_grch38.bed hg38_chm1.fa
