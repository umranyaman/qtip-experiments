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

REF_FASTA_HG=hg19.fa
REF_FASTA_MM=mm10.fa
REF_FASTA_ZM=zm_AGPv3.fa

ifeq (,$(TS_HOME))
$(error "Set TS_HOME, where ts repo root is located")
endif

ifeq (,$(TS_REFS))
$(error "Set TS_REFS, where $(REF_FASTA) is located")
endif

FA_HG=$(TS_REFS)/$(REF_FASTA_HG)
FA_MM=$(TS_REFS)/$(REF_FASTA_MM)
FA_ZM=$(TS_REFS)/$(REF_FASTA_ZM)

ifeq ($(wildcard $(FA_HG)),) 
$(error "$(FA_HG) doesn't exist")
endif
ifeq ($(wildcard $(FA_MM)),) 
$(error "$(FA_MM) doesn't exist")
endif
ifeq ($(wildcard $(FA_ZM)),) 
$(error "$(FA_ZM) doesn't exist")
endif

.PHONY: all
all: reads

include ../mason.mk

.PHONY: clean
clean:
	rm -f *.fq

.PHONY: reads
reads: r0_mason_ill_hg_100_100k.fq.gz \
       r0_mason_ill_hg_250_100k.fq.gz \
       r0_mason_ill_mm_100_100k.fq.gz \
       r0_mason_ill_mm_250_100k.fq.gz \
       r1_mason_ill_hg_100_100k.fq.gz \
       r1_mason_ill_hg_250_100k.fq.gz \
       r1_mason_ill_mm_100_100k.fq.gz \
       r1_mason_ill_mm_250_100k.fq.gz

$(eval $(call mason_ill_unp_reads,mason_ill_hg_100_100k,$(FA_HG),100,100000,7724))
$(eval $(call mason_ill_unp_reads,mason_ill_hg_250_100k,$(FA_HG),250,100000,7724))

$(eval $(call mason_ill_unp_reads,mason_ill_mm_100_100k,$(FA_MM),100,100000,7724))
$(eval $(call mason_ill_unp_reads,mason_ill_mm_250_100k,$(FA_MM),250,100000,7724))

$(eval $(call mason_ill_unp_reads,mason_ill_zm_100_100k,$(FA_ZM),100,100000,7724))
$(eval $(call mason_ill_unp_reads,mason_ill_zm_250_100k,$(FA_ZM),250,100000,7724))

$(eval $(call mason_ill_pair_reads,mason_ill_hg_100_100k,$(FA_HG),100,100000,300,100,7724))
$(eval $(call mason_ill_pair_reads,mason_ill_hg_250_100k,$(FA_HG),250,100000,750,250,7724))

$(eval $(call mason_ill_pair_reads,mason_ill_mm_100_100k,$(FA_MM),100,100000,300,100,7724))
$(eval $(call mason_ill_pair_reads,mason_ill_mm_250_100k,$(FA_MM),250,100000,750,250,7724))

$(eval $(call mason_ill_pair_reads,mason_ill_zm_100_100k,$(FA_ZM),100,100000,300,100,7724))
$(eval $(call mason_ill_pair_reads,mason_ill_zm_250_100k,$(FA_ZM),250,100000,750,250,7724))
