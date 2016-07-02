#!/bin/bash

curl http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz \
    gzip -dc | grep -av '^>' | tee >(tr -cd ACGTN | wc -c | tee mm10_upper.txt | xargs echo "mm10 upper") | \
                                     tr -cd acgtn | wc -c | tee mm10_lower.txt | xargs echo "mm10 lower"

curl http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz \
    gzip -dc | grep -av '^>' | tee >(tr -cd ACGTN | wc -c | tee hg19_upper.txt | xargs echo "hg19 upper") | \
                                     tr -cd acgtn | wc -c | tee hg19_lower.txt | xargs echo "hg19 lower"
