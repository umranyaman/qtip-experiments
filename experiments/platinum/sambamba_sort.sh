#!/bin/sh

NM=final

if [ ! -f "ERR194147.sam/${NM}.bam" ] ; then
    sambamba view -S -f bam ERR194147.sam/${NM}.sam > ERR194147.sam/${NM}.bam
fi
if [ ! -f "ERR194147.sam/${NM}.sorted.bam" ] ; then
    mkdir -p sambamba_temp
    sambamba sort --tmpdir=sambamba_temp -p -m 500G -t 56 ERR194147.sam/${NM}.bam ERR194147.sam/${NM}.sorted.bam
fi
