#!/bin/bash -l
#SBATCH
#SBATCH --job-name=VCFROC
#SBATCH --output=.VCFROC.out
#SBATCH --error=.VCFROC.err
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --partition=shared
#SBATCH --time=4:00:00

# Generate .roc files

VR="$QTIP_EXPERIMENTS_HOME/software/vcflib/vcflib-git/bin/vcfroc"
REFDIR="$QTIP_EXPERIMENTS_HOME/experiments/refs"
NM=ERR194147
SAMP=NA12878

for MINMAPQ in 00 01 02 03 04 05 06 07 08 09 10 11 12 15 20 30 d s u ; do
cat >.VcfRoc.${MINMAPQ}.sh <<EOF
#!/bin/bash -l
#SBATCH
#SBATCH --job-name=VcfRoc
#SBATCH --output=.VcfRoc.${MINMAPQ}.out
#SBATCH --error=.VcfRoc.${MINMAPQ}.err
#SBATCH --nodes=1
#SBATCH --mem=12G
#SBATCH --partition=shared
#SBATCH --time=1:00:00

for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ; do
    for BAM in input final ; do
        FN="${NM}_\${BAM}_\${CHR}_${MINMAPQ}"
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
                    gawk -f f.awk \${ROC_FN} > \${F_FN}
                else
                    echo "Skipping \${ROC_FN}, already present"
                fi
            else
                echo "Skipping \${VCF_FN}, not present"
            fi
        done
    done
done
EOF
echo "sbatch .VcfRoc.${MINMAPQ}.sh"
[ "$1" = "wet" ] && sbatch .VcfRoc.${MINMAPQ}.sh && sleep 1

done
