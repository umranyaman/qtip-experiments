#!/bin/sh

# Generate .roc files

VR="$QTIP_EXPERIMENTS_HOME/software/vcflib/vcflib-git/bin/vcfroc"
REFDIR="$QTIP_EXPERIMENTS_HOME/experiments/refs"
NM=ERR194147
SAMP=NA12878

#for COV in F 50 40 30 ; do
for COV in F ; do
    for MINMAPQ in 00 01 02 03 04 05 06 07 08 09 10 11 12 15 20 30 d s u ; do
        cat >.VcfRoc.${MINMAPQ}.${COV}.sh <<EOF
#!/bin/bash -l
#SBATCH
#SBATCH --job-name=VcfRoc
#SBATCH --output=.VcfRoc.${MINMAPQ}.${COV}.out
#SBATCH --error=.VcfRoc.${MINMAPQ}.${COV}.err
#SBATCH --nodes=1
#SBATCH --mem=12G
#SBATCH --partition=shared
#SBATCH --time=4:00:00

CHR=W
    for BAM in input final ; do
        FN="${NM}_\${BAM}_\${CHR}_${MINMAPQ}_${COV}"
        PREF="${NM}.sam/\${FN}."
        #for FILT in rmsk cr ; do
        for FILT in cr ; do
            VR_COMMON="${VR} --truth-vcf ${SAMP}.\${CHR}.\${FILT}_filt.vcf --reference ${REFDIR}/hg38.fa"
            ROC_FN="\${PREF}\${FILT}_filt.roc"
            VCF_FN="\${PREF}\${FILT}_filt.vcf"
            F_FN="\${PREF}\${FILT}_filt.f"
            if [ -f "\${VCF_FN}" ] ; then
                if [ ! -f "\${ROC_FN}" ] ; then
                    \${VR_COMMON} \${VCF_FN} > \${ROC_FN}
                else
                    echo "Skipping \${ROC_FN}, already present"
                fi
                if [ ! -f "\${F_FN}" ] ; then
                    gawk -f f.awk \${ROC_FN} > \${F_FN}
                else
                    echo "Skipping \${F_FN}, already present"
                fi
            else
                echo "Skipping \${VCF_FN}, not present"
            fi
        done
    done
EOF
        echo "sbatch .VcfRoc.${MINMAPQ}.${COV}.sh"
        [ "$1" = "wet" ] && sbatch .VcfRoc.${MINMAPQ}.${COV}.sh && sleep 1
    done
done
