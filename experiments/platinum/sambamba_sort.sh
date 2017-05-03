#!/bin/bash -l
#SBATCH
#SBATCH --nodes=1
#SBATCH --mem=110G
#SBATCH --partition=parallel
#SBATCH --time=4:00:00

SAMP=ERR194147
NTHREADS=24
MEM="100G"

for NM in input final ; do
    if [ ! -f "${SAMP}.sam/${NM}.bam" ] ; then
	sambamba view -S -f bam "${SAMP}.sam/${NM}.sam" > "${SAMP}.sam/${NM}.bam"
    fi
    if [ ! -f "${SAMP}.sam/${NM}.sorted.bam" ] ; then
	TEMP="sambamba_temp_${SAMP}_${NM}"
	mkdir -p "${TEMP}"
	sambamba sort --tmpdir="${TEMP}" -p -m ${MEM} -t ${NTHREADS} "${SAMP}.sam/${NM}.bam" "${SAMP}.sam/${NM}.sorted.bam"
    fi
done

