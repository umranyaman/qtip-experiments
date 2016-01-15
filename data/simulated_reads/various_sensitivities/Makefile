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
all: outs

BOWTIE2_HOME=$(TS_HOME)/software/bowtie2
BOWTIE2=$(BOWTIE2_HOME)/bowtie2
BT2_ARGS=

$(BOWTIE2):
	$(MAKE) -C $(BOWTIE2_HOME)

include ../mason.mk
include ../bowtie2ts.mk

.PHONY: clean
clean:
	rm -f *.fq.gz
	rm -rf *.out

.PHONY: outs
outs: outs_unpaired outs_paired

.PHONY: outs_unpaired
outs_unpaired: \
    r0_bt2vf_mason_ill_100_1M.out \
    r0_bt2vf_mason_ill_250_1M.out \
    r0_bt2f_mason_ill_100_1M.out \
    r0_bt2f_mason_ill_250_1M.out \
    r0_bt2s_mason_ill_100_1M.out \
    r0_bt2s_mason_ill_250_1M.out \
    r0_bt2vs_mason_ill_100_1M.out \
    r0_bt2vs_mason_ill_250_1M.out \
    r0_bt2vfl_mason_ill_100_1M.out \
    r0_bt2vfl_mason_ill_250_1M.out \
    r0_bt2fl_mason_ill_100_1M.out \
    r0_bt2fl_mason_ill_250_1M.out \
    r0_bt2sl_mason_ill_100_1M.out \
    r0_bt2sl_mason_ill_250_1M.out \
    r0_bt2vsl_mason_ill_100_1M.out \
    r0_bt2vsl_mason_ill_250_1M.out

.PHONY: outs_paired
outs_paired: \
    r12_bt2vf100_mason_ill_100_1M.out \
    r12_bt2vf250_mason_ill_250_1M.out \
    r12_bt2f100_mason_ill_100_1M.out \
    r12_bt2f250_mason_ill_250_1M.out \
    r12_bt2s100_mason_ill_100_1M.out \
    r12_bt2s250_mason_ill_250_1M.out \
    r12_bt2vs100_mason_ill_100_1M.out \
    r12_bt2vs250_mason_ill_250_1M.out \
    r12_bt2vfl100_mason_ill_100_1M.out \
    r12_bt2vfl250_mason_ill_250_1M.out \
    r12_bt2fl100_mason_ill_100_1M.out \
    r12_bt2fl250_mason_ill_250_1M.out \
    r12_bt2sl100_mason_ill_100_1M.out \
    r12_bt2sl250_mason_ill_250_1M.out \
    r12_bt2vsl100_mason_ill_100_1M.out \
    r12_bt2vsl250_mason_ill_250_1M.out

.PHONY: reads
reads: \
    r0_mason_ill_100_1M.fq.gz \
    r0_mason_ill_250_1M.fq.gz \
    r1_mason_ill_100_1M.fq.gz \
    r1_mason_ill_250_1M.fq.gz

SUBSAMPLING_ARGS=--subsampling-series 0.01,0.03,0.05,0.1 --trials=5 --predict-for-training --keep-intermediates --assess-accuracy

#
# Make bowtie2 rule.  Use --sensitive, hg19, and prefix rules with "bt2s"
#
$(eval $(call bt2ts,bt2vf,1000000,1000000,$(SUBSAMPLING_ARGS),--very-fast,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2f,1000000,1000000,$(SUBSAMPLING_ARGS),--fast,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2s,1000000,1000000,$(SUBSAMPLING_ARGS),--sensitive,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2vs,1000000,1000000,$(SUBSAMPLING_ARGS),--very-sensitive,$(HG_FASTA),$(HG_INDEX)))

$(eval $(call bt2ts,bt2vfl,1000000,1000000,$(SUBSAMPLING_ARGS),--very-fast-local,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2fl,1000000,1000000,$(SUBSAMPLING_ARGS),--fast-local,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2sl,1000000,1000000,$(SUBSAMPLING_ARGS),--sensitive-local,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2vsl,1000000,1000000,$(SUBSAMPLING_ARGS),--very-sensitive-local,$(HG_FASTA),$(HG_INDEX)))

$(eval $(call bt2ts,bt2vf100,1000000,1000000,$(SUBSAMPLING_ARGS),-I 200 -X 400 --very-fast,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2vf250,1000000,1000000,$(SUBSAMPLING_ARGS),-I 500 -X 1000 --very-fast,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2f100,1000000,1000000,$(SUBSAMPLING_ARGS),-I 200 -X 400 --fast,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2f250,1000000,1000000,$(SUBSAMPLING_ARGS),-I 500 -X 1000 --fast,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2s100,1000000,1000000,$(SUBSAMPLING_ARGS),-I 200 -X 400 --sensitive,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2s250,1000000,1000000,$(SUBSAMPLING_ARGS),-I 500 -X 1000 --sensitive,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2vs100,1000000,1000000,$(SUBSAMPLING_ARGS),-I 200 -X 400 --very-sensitive,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2vs250,1000000,1000000,$(SUBSAMPLING_ARGS),-I 500 -X 1000 --very-sensitive,$(HG_FASTA),$(HG_INDEX)))

$(eval $(call bt2ts,bt2vfl100,1000000,1000000,$(SUBSAMPLING_ARGS),-I 200 -X 400 --very-fast-local,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2vfl250,1000000,1000000,$(SUBSAMPLING_ARGS),-I 500 -X 1000 --very-fast-local,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2fl100,1000000,1000000,$(SUBSAMPLING_ARGS),-I 200 -X 400 --fast-local,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2fl250,1000000,1000000,$(SUBSAMPLING_ARGS),-I 500 -X 1000 --fast-local,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2sl100,1000000,1000000,$(SUBSAMPLING_ARGS),-I 200 -X 400 --sensitive-local,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2sl250,1000000,1000000,$(SUBSAMPLING_ARGS),-I 500 -X 1000 --sensitive-local,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2vsl100,1000000,1000000,$(SUBSAMPLING_ARGS),-I 200 -X 400 --very-sensitive-local,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2vsl250,1000000,1000000,$(SUBSAMPLING_ARGS),-I 500 -X 1000 --very-sensitive-local,$(HG_FASTA),$(HG_INDEX)))

$(eval $(call mason_ill_unp_reads,ill_100_1M,$(HG_FASTA),100,1000000,7724))
$(eval $(call mason_ill_unp_reads,ill_250_1M,$(HG_FASTA),250,1000000,7724))

$(eval $(call mason_ill_pair_reads,ill_100_1M,$(HG_FASTA),100,1000000,300,100,7724))
$(eval $(call mason_ill_pair_reads,ill_250_1M,$(HG_FASTA),250,1000000,750,250,7724))
