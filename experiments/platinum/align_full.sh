#!/bin/bash

ALIGNER_CPUS=$1
[ -z "${ALIGNER_CPUS}" ] && ALIGNER_CPUS=24

# 1: NAME
# 2: FASTQ1
# 3: FASTQ2
make_job() {
    SCR_FN=".AlignFull.${NM}.sh"
    cat >${SCR_FN} << EOF
#!/bin/bash -l
#SBATCH
#SBATCH --job-name=AlignFull.${NM}
#SBATCH --output=.AlignFull.${NM}.out
#SBATCH --error=.AlignFull.${NM}.err
#SBATCH --nodes=1
#SBATCH --mem=12G
#SBATCH --partition=shared
#SBATCH --cpus-per-task=${ALIGNER_CPUS}
#SBATCH --time=48:00:00

ODIR="${1}.sam"
if [ ! -f "${ODIR}/final.sam" ] ; then
    TEMP="${NM}.temp"
    rm -rf ${TEMP}
    mkdir -p ${TEMP}
    ${QTIP_HOME}/src/qtip \
        --ref ${QTIP_EXPERIMENTS_HOME}/experiments/refs/hg38.fa \
        --m1 ${2} --m2 ${3} \
        --index ${QTIP_EXPERIMENTS_HOME}/experiments/refs/hg38.fa \
        --bt2-exe ${QTIP_HOME}/software/bowtie2/bowtie2 \
        --keep-intermediates \
        --output-directory ${ODIR} \
        --write-orig-mapq \
        --write-precise-mapq \
        --temp-directory ${TEMP} \
        -- -I 0 -X 550 -t -p${ALIGNER_CPUS} --reorder
fi
EOF
    echo "sbatch ${SCR_FN}"
    if [ "${4}" = "wet" ] ; then
        sbatch ${SCR_FN}
    fi
}

for samp in 6 7 ; do
    NM="ERR19414${samp}"
    make_job ${NM} "${NM}_1.fastq.gz" "${NM}_2.fastq.gz" ${1}
done

for samp in 36 37 38 39 40 41 ; do
    NM="SRR6426${samp}"
    make_job ${NM} "${NM}_1.fastq.gz" "${NM}_2.fastq.gz" ${1}
done
