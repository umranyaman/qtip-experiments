#!/bin/bash -l
#SBATCH
#SBATCH --nodes=1
#SBATCH --mem=110G
#SBATCH --partition=parallel
#SBATCH --time=12:00:00

NM=final
SAMP=ERR194147
NTHREADS=24
MEM="100G"

if [ ! -f "${SAMP}.sam/${NM}.bam" ] ; then
    sambamba view -S -f bam "${SAMP}.sam/${NM}.sam" > "${SAMP}.sam/${NM}.bam"
fi
if [ ! -f "${SAMP}.sam/${NM}.sorted.bam" ] ; then
    mkdir -p "sambamba_temp_${SAMP}"
    sambamba sort --tmpdir="sambamba_temp_${SAMP}" -p -m ${MEM} -t ${NTHREADS} "${SAMP}.sam/${NM}.bam" "${SAMP}.sam/${NM}.sorted.bam"
fi
