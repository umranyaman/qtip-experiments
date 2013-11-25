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
include ../wgsim.mk
include ../art.mk
include ../bowtie2.mk

.PHONY: clean
clean:
	rm -f *.fq

.PHONY: alignments
alignments: r0_bt2s_mason_ill_hg_100_100k.sam.gz \
            r0_bt2s_mason_ill_hg_250_100k.sam.gz \
            r0_bt2s_wgsim_ill_hg_100_100k.sam.gz \
            r0_bt2s_wgsim_ill_hg_250_100k.sam.gz \
            r0_bt2s_art_ill_hg_100_100k.sam.gz \
            r0_bt2s_art_ill_hg_250_100k.sam.gz \
            r12_bt2s_mason_ill_hg_100_100k.sam.gz \
            r12_bt2s_mason_ill_hg_250_100k.sam.gz \
            r12_bt2s_wgsim_ill_hg_100_100k.sam.gz \
            r12_bt2s_wgsim_ill_hg_250_100k.sam.gz \
            r12_bt2s_art_ill_hg_100_100k.sam.gz \
            r12_bt2s_art_ill_hg_250_100k.sam.gz

.PHONY: reads
reads: r0_mason_ill_hg_100_100k.fq.gz \
       r0_mason_ill_hg_250_100k.fq.gz \
       r0_wgsim_ill_hg_100_100k.fq.gz \
       r0_wgsim_ill_hg_250_100k.fq.gz \
       r0_art_ill_hg_100_100k.fq.gz \
       r0_art_ill_hg_250_100k.fq.gz \
       r1_mason_ill_hg_100_100k.fq.gz \
       r1_mason_ill_hg_250_100k.fq.gz \
       r1_wgsim_ill_hg_100_100k.fq.gz \
       r1_wgsim_ill_hg_250_100k.fq.gz \
       r1_art_ill_hg_100_100k.fq.gz \
       r1_art_ill_hg_250_100k.fq.gz

#
# Make bowtie2 rule.  Use --sensitive, hg19, and prefix rules with "bt2s"
#
$(eval $(call bt2,bt2s,--sensitive,$(HG_INDEX)))

#
# Make the read rules
#
$(eval $(call mason_ill_unp_reads,ill_hg_100_100k,$(HG_FASTA),100,100000,7724))
$(eval $(call mason_ill_unp_reads,ill_hg_250_100k,$(HG_FASTA),250,100000,7724))
$(eval $(call wgsim_ill_unp_reads,ill_hg_100_100k,$(HG_FASTA),100,100000,7724))
$(eval $(call wgsim_ill_unp_reads,ill_hg_250_100k,$(HG_FASTA),250,100000,7724))

$(eval $(call mason_ill_pair_reads,ill_hg_100_100k,$(HG_FASTA),100,100000,300,25,7724))
$(eval $(call mason_ill_pair_reads,ill_hg_250_100k,$(HG_FASTA),250,100000,750,75,7724))
$(eval $(call wgsim_ill_pair_reads,ill_hg_100_100k,$(HG_FASTA),100,100000,300,25,7724))
$(eval $(call wgsim_ill_pair_reads,ill_hg_250_100k,$(HG_FASTA),250,100000,750,75,7724))

# For Art, we have to convert desired # reads to a fold coverage.
# Note: for paired-end, we're specifying *fragment* coverage.

# >>> 100000 * 100 / 3.1e9
# 0.0032258064516129032
$(eval $(call art_ill_unp_reads,ill_hg_100_100k,$(HG_FASTA),100,0.0033,7724))
# >>> 100000 * 250 / 3.1e9
# 0.008064516129032258
$(eval $(call art_ill_unp_reads,ill_hg_250_100k,$(HG_FASTA),250,0.0081,7724))

# >>> 100000 * 300 / 3.1e9
# 0.00967741935483871
$(eval $(call art_ill_pair_reads,ill_hg_100_100k,$(HG_FASTA),100,0.01,300,25,7724))
# >>> 100000 * 750 / 3.1e9
# 0.024193548387096774
$(eval $(call art_ill_pair_reads,ill_hg_250_100k,$(HG_FASTA),250,0.025,750,75,7724))
