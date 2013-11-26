#!/bin/sh

# Generate simulated FASTQ files for Illumina-like reads, both
# unpaired and paired-end.  We generate several files, each for a
# different read length from 50 up to 500 nt.  Uses the Mason
# simulator.

# Prerequisites:
# 1. TS_HOME (root of repo) must be set
# 2. TS_REFS (containing hg19.fa) must be set
# 3. TS_INDEXES (containing hg19.fa.*.bt2) must be set
# 4. $TS_HOME/software/mason/mason must have been built

HG_FASTA_NAME=hg19.fa
HG_INDEX_BASE=$(HG_FASTA_NAME)

ifeq (,$(TS_HOME))
$(error "Set TS_HOME, where ts repo root is located")
endif
ifeq (,$(TS_REFS))
$(error "Set TS_REFS, where $(REF_FASTA) is located")
endif
ifeq (,$(TS_INDEXES))
$(error "Set TS_INDEXES, where ts repo root is located")
endif

HG_FASTA=$(TS_REFS)/$(HG_FASTA_NAME)
HG_INDEX=$(TS_INDEXES)/$(HG_INDEX_BASE)

ifeq ($(wildcard $(HG_FASTA)),) 
$(error "$(HG_FASTA) doesn't exist")
endif
ifeq ($(wildcard $(HG_INDEX).1.bt2),) 
$(error "$(HG_INDEX).1.bt2 doesn't exist")
endif

.PHONY: all
all: reads

include ../mason.mk
include ../bowtie2.mk

.PHONY: clean
clean:
	rm -f *.fq

.PHONY: alignments
alignments: \
    r0_bt2s_mason_ill_50_100k.sam.gz  \
    r0_bt2s_mason_ill_100_100k.sam.gz \
    r0_bt2s_mason_ill_150_100k.sam.gz \
    r0_bt2s_mason_ill_250_100k.sam.gz \
    r0_bt2s_mason_ill_500_100k.sam.gz \
    r1_bt2s_mason_ill_50_100k.sam.gz  \
    r1_bt2s_mason_ill_100_100k.sam.gz \
    r1_bt2s_mason_ill_150_100k.sam.gz \
    r1_bt2s_mason_ill_250_100k.sam.gz \
    r1_bt2s_mason_ill_500_100k.sam.gz

.PHONY: reads
reads: \
    r0_mason_ill_50_100k.fq.gz  \
    r0_mason_ill_100_100k.fq.gz \
    r0_mason_ill_150_100k.fq.gz \
    r0_mason_ill_250_100k.fq.gz \
    r0_mason_ill_500_100k.fq.gz \
    r1_mason_ill_50_100k.fq.gz  \
    r1_mason_ill_100_100k.fq.gz \
    r1_mason_ill_150_100k.fq.gz \
    r1_mason_ill_250_100k.fq.gz \
    r1_mason_ill_500_100k.fq.gz

#
# Make bowtie2 rule.  Use --sensitive, hg19, and prefix rules with "bt2s"
#
$(eval $(call bt2,bt2s,--sensitive,$(HG_INDEX)))

$(eval $(call mason_ill_unp_reads,ill_50_100k,$(HG_FASTA),50,100000,7724))
$(eval $(call mason_ill_unp_reads,ill_100_100k,$(HG_FASTA),100,100000,7724))
$(eval $(call mason_ill_unp_reads,ill_150_100k,$(HG_FASTA),150,100000,7724))
$(eval $(call mason_ill_unp_reads,ill_250_100k,$(HG_FASTA),250,100000,7724))
$(eval $(call mason_ill_unp_reads,ill_500_100k,$(HG_FASTA),500,100000,7724))

$(eval $(call mason_ill_pair_reads,ill_50_100k,$(HG_FASTA),50,100000,150,50,7724))
$(eval $(call mason_ill_pair_reads,ill_100_100k,$(HG_FASTA),100,100000,300,100,7724))
$(eval $(call mason_ill_pair_reads,ill_150_100k,$(HG_FASTA),150,100000,450,150,7724))
$(eval $(call mason_ill_pair_reads,ill_250_100k,$(HG_FASTA),250,100000,750,250,7724))
$(eval $(call mason_ill_pair_reads,ill_500_100k,$(HG_FASTA),500,100000,1500,500,7724))
