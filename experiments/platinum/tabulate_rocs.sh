#!/bin/sh

NM=ERR194147
awk -v OFS='\t' '{if($1 == "threshold") { if(head == 0) { print $0,"filename" ; head=1 } } else {print $0, FILENAME}}' ${NM}.sam/${NM}_*_*_*.cr_filt.roc > ${NM}_rocs.cr_filt.tsv
