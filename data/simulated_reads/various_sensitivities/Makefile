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
all: tss

BOWTIE2_HOME=$(TS_HOME)/software/bowtie2
BOWTIE2=$(BOWTIE2_HOME)/bowtie2
BT2_ARGS=

$(BOWTIE2):
	$(MAKE) -C $(BOWTIE2_HOME)

include ../mason.mk
include ../bowtie2ts.mk
include ../predict.mk

.PHONY: clean
clean:
	rm -f *.fq.gz
	rm -rf *.ts *.pred

.PHONY: preds
preds: preds_unpaired preds_paired

.PHONY: preds_unpaired
preds_unpaired: \
    r0_bt2vf_mason_ill_100_100k.pred \
    r0_bt2vf_mason_ill_250_100k.pred \
    r0_bt2f_mason_ill_100_100k.pred \
    r0_bt2f_mason_ill_250_100k.pred \
    r0_bt2s_mason_ill_100_100k.pred \
    r0_bt2s_mason_ill_250_100k.pred \
    r0_bt2vs_mason_ill_100_100k.pred \
    r0_bt2vs_mason_ill_250_100k.pred \
    r0_bt2vfl_mason_ill_100_100k.pred \
    r0_bt2vfl_mason_ill_250_100k.pred \
    r0_bt2fl_mason_ill_100_100k.pred \
    r0_bt2fl_mason_ill_250_100k.pred \
    r0_bt2sl_mason_ill_100_100k.pred \
    r0_bt2sl_mason_ill_250_100k.pred \
    r0_bt2vsl_mason_ill_100_100k.pred \
    r0_bt2vsl_mason_ill_250_100k.pred

.PHONY: preds_paired
preds_paired: \
    r12_bt2vf100_mason_ill_100_100k.pred \
    r12_bt2vf250_mason_ill_250_100k.pred \
    r12_bt2f100_mason_ill_100_100k.pred \
    r12_bt2f250_mason_ill_250_100k.pred \
    r12_bt2s100_mason_ill_100_100k.pred \
    r12_bt2s250_mason_ill_250_100k.pred \
    r12_bt2vs100_mason_ill_100_100k.pred \
    r12_bt2vs250_mason_ill_250_100k.pred \
    r12_bt2vfl100_mason_ill_100_100k.pred \
    r12_bt2vfl250_mason_ill_250_100k.pred \
    r12_bt2fl100_mason_ill_100_100k.pred \
    r12_bt2fl250_mason_ill_250_100k.pred \
    r12_bt2sl100_mason_ill_100_100k.pred \
    r12_bt2sl250_mason_ill_250_100k.pred \
    r12_bt2vsl100_mason_ill_100_100k.pred \
    r12_bt2vsl250_mason_ill_250_100k.pred

.PHONY: tss
tss: tss_unpaired tss_paired

.PHONY: tss_unpaired
tss_unpaired: \
    r0_bt2vf_mason_ill_100_100k.ts \
    r0_bt2vf_mason_ill_250_100k.ts \
    r0_bt2f_mason_ill_100_100k.ts \
    r0_bt2f_mason_ill_250_100k.ts \
    r0_bt2s_mason_ill_100_100k.ts \
    r0_bt2s_mason_ill_250_100k.ts \
    r0_bt2vs_mason_ill_100_100k.ts \
    r0_bt2vs_mason_ill_250_100k.ts \
    r0_bt2vfl_mason_ill_100_100k.ts \
    r0_bt2vfl_mason_ill_250_100k.ts \
    r0_bt2fl_mason_ill_100_100k.ts \
    r0_bt2fl_mason_ill_250_100k.ts \
    r0_bt2sl_mason_ill_100_100k.ts \
    r0_bt2sl_mason_ill_250_100k.ts \
    r0_bt2vsl_mason_ill_100_100k.ts \
    r0_bt2vsl_mason_ill_250_100k.ts

.PHONY: tss_paired
tss_paired: \
    r12_bt2vf100_mason_ill_100_100k.ts \
    r12_bt2vf250_mason_ill_250_100k.ts \
    r12_bt2f100_mason_ill_100_100k.ts \
    r12_bt2f250_mason_ill_250_100k.ts \
    r12_bt2s100_mason_ill_100_100k.ts \
    r12_bt2s250_mason_ill_250_100k.ts \
    r12_bt2vs100_mason_ill_100_100k.ts \
    r12_bt2vs250_mason_ill_250_100k.ts \
    r12_bt2vfl100_mason_ill_100_100k.ts \
    r12_bt2vfl250_mason_ill_250_100k.ts \
    r12_bt2fl100_mason_ill_100_100k.ts \
    r12_bt2fl250_mason_ill_250_100k.ts \
    r12_bt2sl100_mason_ill_100_100k.ts \
    r12_bt2sl250_mason_ill_250_100k.ts \
    r12_bt2vsl100_mason_ill_100_100k.ts \
    r12_bt2vsl250_mason_ill_250_100k.ts

.PHONY: reads
reads: \
    r0_mason_ill_100_100k.fq.gz \
    r0_mason_ill_250_100k.fq.gz \
    r1_mason_ill_100_100k.fq.gz \
    r1_mason_ill_250_100k.fq.gz

#
# Make bowtie2 rule.  Use --sensitive, hg19, and prefix rules with "bt2s"
#
$(eval $(call bt2ts,bt2vf,100000,100000,--low-score-bias=7,--very-fast,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2f,100000,100000,--low-score-bias=7,--fast,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2s,100000,100000,--low-score-bias=7,--sensitive,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2vs,100000,100000,--low-score-bias=7,--very-sensitive,$(HG_FASTA),$(HG_INDEX)))

$(eval $(call bt2ts,bt2vfl,100000,100000,,--very-fast-local,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2fl,100000,100000,,--fast-local,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2sl,100000,100000,,--sensitive-local,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2vsl,100000,100000,,--very-sensitive-local,$(HG_FASTA),$(HG_INDEX)))

$(eval $(call bt2ts,bt2vf100,100000,100000,,-I 200 -X 400 --very-fast,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2vf250,100000,100000,,-I 500 -X 1000 --very-fast,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2f100,100000,100000,,-I 200 -X 400 --fast,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2f250,100000,100000,,-I 500 -X 1000 --fast,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2s100,100000,100000,,-I 200 -X 400 --sensitive,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2s250,100000,100000,,-I 500 -X 1000 --sensitive,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2vs100,100000,100000,,-I 200 -X 400 --very-sensitive,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2vs250,100000,100000,,-I 500 -X 1000 --very-sensitive,$(HG_FASTA),$(HG_INDEX)))

$(eval $(call bt2ts,bt2vfl100,100000,100000,,-I 200 -X 400 --very-fast-local,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2vfl250,100000,100000,,-I 500 -X 1000 --very-fast-local,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2fl100,100000,100000,,-I 200 -X 400 --fast-local,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2fl250,100000,100000,,-I 500 -X 1000 --fast-local,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2sl100,100000,100000,,-I 200 -X 400 --sensitive-local,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2sl250,100000,100000,,-I 500 -X 1000 --sensitive-local,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2vsl100,100000,100000,,-I 200 -X 400 --very-sensitive-local,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2vsl250,100000,100000,,-I 500 -X 1000 --very-sensitive-local,$(HG_FASTA),$(HG_INDEX)))

$(eval $(call mason_ill_unp_reads,ill_100_100k,$(HG_FASTA),100,100000,7724))
$(eval $(call mason_ill_unp_reads,ill_250_100k,$(HG_FASTA),250,100000,7724))

$(eval $(call mason_ill_pair_reads,ill_100_100k,$(HG_FASTA),100,100000,300,100,7724))
$(eval $(call mason_ill_pair_reads,ill_250_100k,$(HG_FASTA),250,100000,750,250,7724))
