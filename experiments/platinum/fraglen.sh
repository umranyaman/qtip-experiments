#!/bin/bash

ALIGNER_CPUS=$1
[ -z "${ALIGNER_CPUS}" ] && ALIGNER_CPUS=24

# 1: NAME
# 2: FASTQ1
# 3: FASTQ2
make_job() {
    SCR_FN=".FraglenFull.${NM}.sh"
    cat >${SCR_FN} << EOF
#!/bin/bash -l
#SBATCH
#SBATCH --job-name=FraglenFull.${NM}
#SBATCH --output=.FraglenFull.${NM}.out
#SBATCH --error=.FraglenFull.${NM}.err
#SBATCH --nodes=1
#SBATCH --mem=12G
#SBATCH --partition=shared
#SBATCH --cpus-per-task=${ALIGNER_CPUS}
#SBATCH --time=1:00:00

OUTFN="${1}.fraglen.sam"
if [ ! -f "\${OUTFN}" ] ; then
    ${QTIP_HOME}/software/bowtie2/bowtie2 \
        -x ${QTIP_EXPERIMENTS_HOME}/experiments/refs/hg38.fa \
        -1 ${2} -2 ${3} \
        -I 0 -X 2000 \
        -t \
        -s 10000000 -u 10000000 \
        -p ${ALIGNER_CPUS} \
        -S \${OUTFN}

    awk -v FS='\t' '\$9 > 0 && \$1 !~ /^@/ && \$9 < 10000 {h[int(\$9/10)] += 1} END {for(d in h) {print d*10,h[d]}}' \
        \${OUTFN} | sort -r -n -k1,1 > ${1}.fraglen
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
