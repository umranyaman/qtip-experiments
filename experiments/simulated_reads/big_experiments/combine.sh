#!/bin/bash -l
#SBATCH
#SBATCH --nodes=1
#SBATCH --mem=1G
#SBATCH --partition=shared
#SBATCH --time=4:00:00

if [ ! -f r1_mason_ill_250_50M.fq.gz ] ; then
    gzip -dc r1_mason_ill_250_?_50M.fq.gz | head -n 200000000 | gzip -c > r1_mason_ill_250_50M.fq.gz
fi
if [ ! -f r2_mason_ill_250_50M.fq.gz ] ; then
    gzip -dc r2_mason_ill_250_?_50M.fq.gz | head -n 200000000 | gzip -c > r2_mason_ill_250_50M.fq.gz
fi
if [ ! -f r1_mason_ill_100_50M.fq.gz ] ; then
    gzip -dc r1_mason_ill_100_?_50M.fq.gz | head -n 200000000 | gzip -c > r1_mason_ill_100_50M.fq.gz
fi
if [ ! -f r2_mason_ill_100_50M.fq.gz ] ; then
    gzip -dc r2_mason_ill_100_?_50M.fq.gz | head -n 200000000 | gzip -c > r2_mason_ill_100_50M.fq.gz
fi
