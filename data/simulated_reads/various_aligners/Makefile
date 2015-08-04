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
ifeq ($(wildcard $(HG_INDEX).amb),) 
$(error "$(HG_INDEX).amb doesn't exist")
endif
ifeq ($(wildcard $(HG_INDEX).mosaik),)
$(error "$(HG_INDEX).mosaik doesn't exist")
endif
ifeq ($(wildcard $(HG_INDEX).snap),)
$(error "$(HG_INDEX).snap doesn't exist")
endif

.PHONY: all
all: tss

BOWTIE2_HOME=$(TS_HOME)/software/bowtie2
BOWTIE2=$(BOWTIE2_HOME)/bowtie2
BT2_ARGS=

BWA_HOME=$(TS_HOME)/software/bwa
BWA=$(BWA_HOME)/bwa
BWA_ARGS=

SNAP_HOME=$(TS_HOME)/software/snap
SNAP=$(SNAP_HOME)/snap/snap-aligner
SNAP_ARGS=

MOSAIK_HOME=$(TS_HOME)/software/mosaik
MOSAIK_ALIGN=$(MOSAIK_HOME)/MosaikAligner
MOSAIK_BUILD=$(MOSAIK_HOME)/MosaikBuild
MOSAIK_ALIGN_ARGS=
MOSAIK_BUILD_ARGS=

$(BOWTIE2):
	$(MAKE) -C $(BOWTIE2_HOME)

$(BWA):
	$(MAKE) -C $(BWA_HOME)

$(SNAP):
	$(MAKE) -C $(BWA_HOME)

$(MOSAIK_ALIGN) $(MOSAIK_BUILD):
	$(MAKE) -C $(MOSAIK_HOME)


include ../mason.mk
include ../bowtie2ts.mk
include ../bwamemts.mk
include ../snapts.mk
include ../mosaikts.mk
include ../predict.mk

.PHONY: clean
clean:
	rm -f *.fq.gz
	rm -f *.mkb
	rm -rf *.ts

.PHONY: preds
preds: preds_bt2 preds_bwamem preds_snap preds_mosaik

.PHONY: preds_bt2
preds_bt2: \
    r0_bt2s_mason_ill_100_100k.pred \
    r0_bt2s_mason_ill_250_100k.pred \
    r12_bt2s100_mason_ill_100_100k.pred \
    r12_bt2s250_mason_ill_250_100k.pred \
    r0_bt2sl_mason_ill_100_100k.pred \
    r0_bt2sl_mason_ill_250_100k.pred \
    r12_bt2sl100_mason_ill_100_100k.pred \
    r12_bt2sl250_mason_ill_250_100k.pred

.PHONY: preds_bwamem
preds_bwamem: \
    r0_bwamem_mason_ill_100_100k.pred \
    r0_bwamem_mason_ill_250_100k.pred \
    r12_bwamem100_mason_ill_100_100k.pred \
    r12_bwamem250_mason_ill_250_100k.pred

.PHONY: preds_mosaik
preds_mosaik: \
    r0_mosaik_mason_ill_100_100k.pred \
    r0_mosaik_mason_ill_250_100k.pred

.PHONY: preds_snap
preds_snap: \
    r0_snap_mason_ill_100_100k.pred \
    r0_snap_mason_ill_250_100k.pred

.PHONY: tss
tss: tss_bt2 tss_bwamem tss_snap tss_mosaik

.PHONY: tss_bt2
tss_bt2: \
    r0_bt2s_mason_ill_100_100k.ts \
    r0_bt2s_mason_ill_250_100k.ts \
    r12_bt2s100_mason_ill_100_100k.ts \
    r12_bt2s250_mason_ill_250_100k.ts \
    r0_bt2sl_mason_ill_100_100k.ts \
    r0_bt2sl_mason_ill_250_100k.ts \
    r12_bt2sl100_mason_ill_100_100k.ts \
    r12_bt2sl250_mason_ill_250_100k.ts

.PHONY: tss_bwamem
tss_bwamem: \
    r0_bwamem_mason_ill_100_100k.ts \
    r0_bwamem_mason_ill_250_100k.ts \
    r12_bwamem100_mason_ill_100_100k.ts \
    r12_bwamem250_mason_ill_250_100k.ts

.PHONY: tss_mosaik
tss_mosaik: \
    r0_mosaik_mason_ill_100_100k.ts \
    r0_mosaik_mason_ill_250_100k.ts

.PHONY: tss_snap
tss_snap: \
    r0_snap_mason_ill_100_100k.ts \
    r0_snap_mason_ill_250_100k.ts

.PHONY: reads
reads: \
    r0_mason_ill_100_100k.fq.gz \
    r0_mason_ill_250_100k.fq.gz \
    r1_mason_ill_100_100k.fq.gz \
    r1_mason_ill_250_100k.fq.gz

#
# Make bowtie2 rule.  Use --sensitive, hg19, and prefix rules with "bt2s"
#
$(eval $(call bt2ts,bt2s,100000,100000,,--sensitive,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2sl,100000,100000,,--sensitive --local,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2s100,100000,100000,,-I 200 -X 400 --sensitive,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2sl100,100000,100000,,-I 200 -X 400 --sensitive --local,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2s250,100000,100000,,-I 500 -X 1000 --sensitive,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2sl250,100000,100000,,-I 500 -X 1000 --sensitive --local,$(HG_FASTA),$(HG_INDEX)))

, := ,

$(eval $(call bwamemts,bwamem,100000,100000,,,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bwamemts,bwamem100,100000,100000,,-I 300$(,)25,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bwamemts,bwamem250,100000,100000,,-I 750$(,)62,$(HG_FASTA),$(HG_INDEX)))

$(eval $(call mosaikts,mosaik,100000,100000,,-a fast,$(HG_FASTA),$(HG_INDEX)))
#$(eval $(call mosaikts,mosaik100,100000,100000,,-I 300$(,)25,$(HG_FASTA),$(HG_INDEX)))
#$(eval $(call mosaikts,mosaik250,100000,100000,,-I 750$(,)62,$(HG_FASTA),$(HG_INDEX)))

$(eval $(call snapts,snap,100000,100000,,,,,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call snapts,snapts100,100000,100000,,,,-s 100 300,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call snapts,snapts250,100000,100000,,,,-s 250 750,$(HG_FASTA),$(HG_INDEX)))

$(eval $(call mason_ill_unp_reads,ill_100_100k,$(HG_FASTA),100,100000,7724))
$(eval $(call mason_ill_unp_reads,ill_250_100k,$(HG_FASTA),250,100000,7724))

$(eval $(call mason_ill_pair_reads,ill_100_100k,$(HG_FASTA),100,100000,300,100,7724))
$(eval $(call mason_ill_pair_reads,ill_250_100k,$(HG_FASTA),250,100000,750,250,7724))
