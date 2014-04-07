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
	rm -rf *.ts

.PHONY: preds
preds: preds_unpaired preds_paired

.PHONY: preds_unpaired
preds_unpaired: \
    r0_bt2s_mason_ill_250_1_100k.pred \
    r0_bt2s_mason_ill_250_2_100k.pred \
    r0_bt2s_mason_ill_250_3_100k.pred \
    r0_bt2s_mason_ill_250_4_100k.pred \
    r0_bt2s_mason_ill_250_5_100k.pred

.PHONY: preds_paired
preds_paired: \
    r12_bt2s100_mason_ill_100_1_100k.pred \
    r12_bt2s100_mason_ill_100_2_100k.pred \
    r12_bt2s100_mason_ill_100_3_100k.pred \
    r12_bt2s100_mason_ill_100_4_100k.pred \
    r12_bt2s100_mason_ill_100_5_100k.pred

.PHONY: tss_unpaired
tss_unpaired: \
    r0_bt2s_mason_ill_250_1_100k.ts \
    r0_bt2s_mason_ill_250_2_100k.ts \
    r0_bt2s_mason_ill_250_3_100k.ts \
    r0_bt2s_mason_ill_250_4_100k.ts \
    r0_bt2s_mason_ill_250_5_100k.ts

.PHONY: tss_paired
tss_paired: \
    r12_bt2s100_mason_ill_100_1_100k.ts \
    r12_bt2s100_mason_ill_100_2_100k.ts \
    r12_bt2s100_mason_ill_100_3_100k.ts \
    r12_bt2s100_mason_ill_100_4_100k.ts \
    r12_bt2s100_mason_ill_100_5_100k.ts

.PHONY: tss
tss: tss_unpaired tss_paired

.PHONY: tss_unpaired
tss_unpaired: \
    r0_bt2s_mason_ill_250_1_100k.ts \
    r0_bt2s_mason_ill_250_2_100k.ts \
    r0_bt2s_mason_ill_250_3_100k.ts \
    r0_bt2s_mason_ill_250_4_100k.ts \
    r0_bt2s_mason_ill_250_5_100k.ts

.PHONY: tss_paired
tss_paired: \
    r12_bt2s100_mason_ill_100_1_100k.ts \
    r12_bt2s100_mason_ill_100_2_100k.ts \
    r12_bt2s100_mason_ill_100_3_100k.ts \
    r12_bt2s100_mason_ill_100_4_100k.ts \
    r12_bt2s100_mason_ill_100_5_100k.ts

# Reads
.PHONY: reads
reads: \
    r0_mason_ill_250_1_100k.fq.gz \
    r0_mason_ill_250_2_100k.fq.gz \
    r0_mason_ill_250_3_100k.fq.gz \
    r0_mason_ill_250_4_100k.fq.gz \
    r0_mason_ill_250_5_100k.fq.gz \
    r1_mason_ill_100_1_100k.fq.gz \
    r1_mason_ill_100_2_100k.fq.gz \
    r1_mason_ill_100_3_100k.fq.gz \
    r1_mason_ill_100_4_100k.fq.gz \
    r1_mason_ill_100_5_100k.fq.gz

#
# Make bowtie2 rule.  Use --sensitive, hg19, and prefix rules with "bt2s"
#
$(eval $(call bt2ts,bt2s,100000,100000,,--sensitive,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2s100,100000,100000,,-I 200 -X 400 --sensitive,$(HG_FASTA),$(HG_INDEX)))

$(eval $(call mason_ill_unp_reads,ill_250_1_100k,$(HG_FASTA),250,100000,7724))
$(eval $(call mason_ill_unp_reads,ill_250_2_100k,$(HG_FASTA),250,100000,1523))
$(eval $(call mason_ill_unp_reads,ill_250_3_100k,$(HG_FASTA),250,100000,1493))
$(eval $(call mason_ill_unp_reads,ill_250_4_100k,$(HG_FASTA),250,100000,9295))
$(eval $(call mason_ill_unp_reads,ill_250_5_100k,$(HG_FASTA),250,100000,6018))

$(eval $(call mason_ill_pair_reads,ill_100_1_100k,$(HG_FASTA),100,100000,300,100,7724))
$(eval $(call mason_ill_pair_reads,ill_100_2_100k,$(HG_FASTA),100,100000,300,100,1523))
$(eval $(call mason_ill_pair_reads,ill_100_3_100k,$(HG_FASTA),100,100000,300,100,1493))
$(eval $(call mason_ill_pair_reads,ill_100_4_100k,$(HG_FASTA),100,100000,300,100,9295))
$(eval $(call mason_ill_pair_reads,ill_100_5_100k,$(HG_FASTA),100,100000,300,100,6018))
