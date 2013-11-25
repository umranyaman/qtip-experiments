#!/bin/sh

# Generate simulated FASTQ files for Illumina-like reads, both
# unpaired and paired-end.  We generate several files, each for a
# different read length from 50 up to 500 nt.  Uses the Mason
# simulator.

# Prerequisites:
# 1. TS_HOME (root of repo) must be set
# 2. TS_REFS (contains hg19/hg19.fa) must be set
# 3. $(TS_REFS)/hg19/hg19.fa must exist
# 4. $TS_HOME/software/mason/mason must have been built

REF_FASTA=hg19.fa

ifeq (,$(TS_HOME))
$(error "Set TS_HOME, where ts repo root is located")
endif

ifeq (,$(TS_REFS))
$(error "Set TS_REFS, where $(REF_FASTA) is located")
endif

FA=$(TS_REFS)/$(REF_FASTA)

ifeq ($(wildcard $(FA)),) 
$(error "$(FA) doesn't exist")
endif

.PHONY: all
all: reads

include ../mason.mk

.PHONY: clean
clean:
	rm -f *.fq

.PHONY: reads
reads: r0_mason_fff_50to500_100k.fq.gz

$(eval $(call mason_fff_unp_reads,mason_fff_50to500_100k,$(FA),275,225,100000,7724))
