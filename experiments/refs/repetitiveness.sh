#!/bin/bash

mkdir -p mm10
cd mm10
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz
tar xvf chromFa.tar.gz
cat *.fa | grep -av '^>' | tee >(tr -cd ACGTN | wc -c | tee ../mm10_upper.txt | xargs echo "mm10 upper") | \
                                 tr -cd acgtn | wc -c | tee ../mm10_lower.txt | xargs echo "mm10 lower"
cd ..
rm -rf mm10

mkdir -p hg19
cd hg19
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
tar xvf chromFa.tar.gz
cat *.fa | grep -av '^>' | tee >(tr -cd ACGTN | wc -c | tee ../hg19_upper.txt | xargs echo "hg19 upper") | \
                                 tr -cd acgtn | wc -c | tee ../hg19_lower.txt | xargs echo "hg19 lower"
cd ..
rm -rf hg19

mkdir -p hg38
cd hg38
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gzip -dc hg38.fa.gz | grep -av '^>' | tee >(tr -cd ACGTN | wc -c | tee ../hg38_upper.txt | xargs echo "hg38 upper") | \
                                            tr -cd acgtn | wc -c | tee ../hg38_lower.txt | xargs echo "hg38 lower"
cd ..
rm -rf hg38
