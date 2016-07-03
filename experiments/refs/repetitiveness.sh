#!/bin/bash

mkdir -p mm10
cd mm10
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz
tar xvf chromfa.tar.gz
cat *.fa | grep -av '^>' | tee >(tr -cd acgtn | wc -c | tee mm10_upper.txt | xargs echo "mm10 upper") | \
                                 tr -cd acgtn | wc -c | tee mm10_lower.txt | xargs echo "mm10 lower"
cd ..
rm -rf mm10

mkdir -p hg19
cd hg19
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
tar xvf chromfa.tar.gz
cat *.fa | grep -av '^>' | tee >(tr -cd acgtn | wc -c | tee mm10_upper.txt | xargs echo "hg19 upper") | \
                                 tr -cd acgtn | wc -c | tee mm10_lower.txt | xargs echo "hg19 lower"
cd ..
rm -rf hg19
