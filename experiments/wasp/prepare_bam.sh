#!/bin/bash -l

#SBATCH
#SBATCH --job-name=PrepareSam
#SBATCH --time=24:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1

samtools sort ../real_data/ERR050082_1.bt2.unp.sam -o ERR050082_1.bt2.unp.sorted.bam
