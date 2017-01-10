#!/bin/bash -l
#SBATCH
#SBATCH --nodes=1
#SBATCH --mem=1G
#SBATCH --partition=shared
#SBATCH --time=3:00:00

cd ERR194147.sam && samtools flagstat final.sorted.bam | tee final.sorted.flagstat
