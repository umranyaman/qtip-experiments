#!/usr/bin/env bash

# Human reads from 1000 Genomes project
# ERR050082 has 42,245,074 paired-end 100 x 100 nt reads
# ERR050083 has 66,391,067 paired-end 100 x 100 nt reads

for d in ERR050082 ERR050083 ; do

    pref=`echo ${d} | sed 's/\(......\).*/\1/'`

    for m in 1 2 ; do
        if [ ! -f "${d}_${m}.fastq" ] ; then
            curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${pref}/${d}/${d}_${m}.fastq.gz | gzip -dc > ${d}_${m}.fastq
        fi
    done
done

for m in 1 2 ; do
    curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/009/SRR1592709/SRR1592709_${m}.fastq.gz | gzip -dc > SRR1592709_${m}.fastq
done

# Zea mays reads