#!/bin/sh

VCFLIB_HOME=$QTIP_EXPERIMENTS_HOME/software/vcflib
VCFISECT=$VCFLIB_HOME/vcflib-git/bin/vcfintersect
NM="NA12878"
PLAT_GET="wget --user=platgene_ro --password= "

[ ! -x "${VCFISECT}" ] && echo "Need vcfintersect binary" && exit 1

# Get gold-standard variant calls
if [ ! -f "${NM}.vcf.gz" ] ; then
    ${PLAT_GET} ftp://platgene_ro@ussd-ftp.illumina.com/2016-1.0/hg38/small_variants/${NM}/${NM}.vcf.gz
fi

# Get info about Platinum confidence
for FN in ConfidentRegions.bed.gz ConfidentRegions.bed.gz.tbi ; do
    if [ ! -f "${FN}" ] ; then
        ${PLAT_GET} ftp://platgene_ro@ussd-ftp.illumina.com/2016-1.0/hg38/small_variants/${FN}
    fi
done

for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ; do
    #[ ! -f "rmsk_${CHR}.bed" ] && echo "Need rmsk file for ${CHR}; see get_low_complexity.sh" && exit 1
    if [ ! -f "cr_${CHR}.bed" ] ; then
        gzip -dc ConfidentRegions.bed.gz | awk "\$1 == \"chr${CHR}\" {print ${CHR},\$2,\$3}" > cr_${CHR}.bed
    fi
    if [ ! -f "${NM}.${CHR}.raw.vcf" ] ; then
        # make chromosome-specific VCF file with normalized chromosome name
        gzip -dc ${NM}.vcf.gz | awk "\$1 ~ /^#/ || \$1 == \"chr${CHR}\"" | sed 's/^chr//' > ${NM}.${CHR}.raw.vcf
    fi
    #if [ ! -f "${NM}.${CHR}.rmsk_filt.vcf" ] ; then
    ## make filtered VCF file with low complexity regions removed
    #${VCFISECT} -v -b rmsk_${CHR}.bed ${NM}.${CHR}.raw.vcf > ${NM}.${CHR}.rmsk_filt.vcf
    #fi
    if [ ! -f "${NM}.${CHR}.cr_filt.vcf" ] ; then
        # make filtered VCF file with Platinum low-confidence regions removed
        ${VCFISECT} -b cr_${CHR}.bed ${NM}.${CHR}.raw.vcf > ${NM}.${CHR}.cr_filt.vcf
    fi
done
