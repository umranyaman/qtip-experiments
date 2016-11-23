#!/bin/sh

#SBATCH
#SBATCH --job-name=Wasp
#SBATCH --time=24:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1

WASP_VER=0.2.1
WASP_DIR="../../software/wasp/WASP-${WASP_VER}"
FIND_SNPS="python ${WASP_DIR}/mapping/find_intersecting_snps.py"

mkdir -p out
$FIND_SNPS --is_sorted --output_dir out --snp_dir snps ERR050082_1.bt2.unp.sorted.bam
