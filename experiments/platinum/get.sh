#!/bin/sh

# Same as the reads used in Heng Li's paper "Toward better
# understanding of artifacts in variant calling from high-coverage
# samples"

for RUN in ERR194146 ERR194147 ; do
    if [ ! -f "${RUN}.fastq.gz" ] ; then
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/${RUN}/${RUN}.fastq.gz
    fi
    if [ ! -f "${RUN}_1.fastq.gz" ] ; then
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/${RUN}/${RUN}_1.fastq.gz
    fi
    if [ ! -f "${RUN}_2.fastq.gz" ] ; then
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/${RUN}/${RUN}_2.fastq.gz
    fi
done

for samp in 36 37 38 39 40 41 ; do
    if [ ! -f "SRR6426${samp}_1.fastq.gz" ] ; then
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR642/SRR6426${samp}/SRR6426${samp}_1.fastq.gz
    fi
    if [ ! -f "SRR6426${samp}_2.fastq.gz" ] ; then
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR642/SRR6426${samp}/SRR6426${samp}_2.fastq.gz
    fi
done
