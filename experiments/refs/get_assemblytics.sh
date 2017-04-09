#!/bin/sh

#
# Make FASTA files that include unaligned/poorly-aligned sequence from CHM1 assembies
# Both are MHAP assemblies
#

get_assembly() {
    if [ ! -f ${1} ] ; then
	curl ${2} | gzip -dc > ${1}
	samtools faidx ${1}
	awk -v OFS='\t' '{print $1,$2}' "${1}.fai" > "${1}.genome"
    fi
}

get_assembly chm1_quiver.fa http://gembox.cbcb.umd.edu/mhap/asm/human.quiver.ctg.fasta.gz
get_assembly chm1_mhap.fa http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/genomes/genbank/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCA_000772585.3_ASM77258v3/GCA_000772585.3_ASM77258v3_genomic.fna.gz

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
        mv user_data $1
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
        # 5. Filters out intervals shorter than 260 bp (for simulation)
        awk '{n += 1; if(n > 1) {if($7 == "Insertion") {print $10}}}' ${BED} | \
	    sed 's/[:-]/,/g' | \
	    awk -v FS=',' -v OFS='\t' '$3-$2 >= 260 {print $1,$2,$3}' > .tmp.bed
        # Add slop on both sides
        bedtools slop -i .tmp.bed -g "${6}.genome" -b ${SLOP} > "${4}"
    fi

    # Now extract sequences from CHM1
    if [ ! -f "${5}" ] ; then
        bedtools getfasta -fi "${6}" -bed "${4}" -fo "${5}"
    fi
}

getbed ${BED_37} 25 50 chm1_grch19.bed hg19_chm1.fa chm1_mhap.fa
getbed ${BED_38} 25 50 chm1_grch38.bed hg38_chm1.fa chm1_quiver.fa

cat hg19.fa hg19_chm1.fa > hg19_with_chm1.fa
cat hg38.fa hg38_chm1.fa > hg38_with_chm1.fa
