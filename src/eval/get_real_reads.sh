#!/usr/bin/env bash

for d in ERR050082 ERR050083 ; do

    pref=`echo ${d} | sed 's/\(......\).*/\1/'`

    for m in 1 2 ; do
        curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${pref}/${d}/${d}_${m}.fastq.gz | gzip -dc > ${d}_${m}.fastq
    done
done
