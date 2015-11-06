#!/usr/bin/env bash

curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR050/ERR050082/ERR050082_1.fastq.gz | gzip -dc | head -n 10000000 > ERR050082_1.excerpt.fastq
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR050/ERR050082/ERR050082_2.fastq.gz | gzip -dc | head -n 10000000 > ERR050082_2.excerpt.fastq

curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR050/ERR050083/ERR050083_1.fastq.gz | gzip -dc | head -n 10000000 > ERR050083_1.excerpt.fastq
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR050/ERR050083/ERR050083_2.fastq.gz | gzip -dc | head -n 10000000 > ERR050083_2.excerpt.fastq
