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

ifeq (,$(TS_HOME))
$(error "Set TS_HOME, where ts repo root is located")
endif

ifeq (,$(TS_REFS))
$(error "Set TS_REFS, where $(REF_FASTA) is located")
endif

FA_HG=$(TS_REFS)/$(REF_FASTA_HG)

ifeq ($(wildcard $(FA_HG)),) 
$(error "$(FA_HG) doesn't exist")
endif

.PHONY: all
all: reads

include ../mason.mk
include ../wgsim.mk
include ../art.mk

.PHONY: clean
clean:
	rm -f *.fq

.PHONY: reads
reads: r0_mason_ill_hg_100_100k.fq.gz \
       r0_mason_ill_hg_250_100k.fq.gz \
       r0_wgsim_ill_hg_100_100k.fq.gz \
       r0_wgsim_ill_hg_250_100k.fq.gz \
       r0_art_ill_hg_100_100k.fq.gz \
       r0_art_ill_hg_250_100k.fq.gz

$(eval $(call mason_ill_unp_reads,ill_hg_100_100k,$(FA_HG),100,100000,7724))
$(eval $(call mason_ill_unp_reads,ill_hg_250_100k,$(FA_HG),250,100000,7724))

$(eval $(call wgsim_ill_unp_reads,ill_hg_100_100k,$(FA_HG),100,100000,7724))
$(eval $(call wgsim_ill_unp_reads,ill_hg_250_100k,$(FA_HG),250,100000,7724))

# >>> 100000 * 100 / 3.1e9
# 0.0032258064516129032
$(eval $(call art_ill_unp_reads,ill_hg_100_100k,$(FA_HG),100,0.0033,7724))
# >>> 100000 * 250 / 3.1e9
# 0.008064516129032258
$(eval $(call art_ill_unp_reads,ill_hg_250_100k,$(FA_HG),250,0.0081,7724))


$(eval $(call mason_ill_pair_reads,ill_hg_100_100k,$(FA_HG),100,100000,300,100,7724))
$(eval $(call mason_ill_pair_reads,ill_hg_250_100k,$(FA_HG),250,100000,750,250,7724))

$(eval $(call wgsim_ill_pair_reads,ill_hg_100_100k,$(FA_HG),100,100000,300,25,7724))
$(eval $(call wgsim_ill_pair_reads,ill_hg_250_100k,$(FA_HG),250,100000,750,75,7724))

# >>> 100000 * 300 / 3.1e9
# 0.00967741935483871
$(eval $(call art_ill_pair_reads,ill_hg_100_100k,$(FA_HG),100,0.01,300,25,7724))
# >>> 100000 * 750 / 3.1e9
# 0.024193548387096774
$(eval $(call art_ill_pair_reads,ill_hg_250_100k,$(FA_HG),250,0.025,750,75,7724))
