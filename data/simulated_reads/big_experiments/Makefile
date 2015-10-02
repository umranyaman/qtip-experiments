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

$(BOWTIE2):
	$(MAKE) -C $(BOWTIE2_HOME)

$(BWA):
	$(MAKE) -C $(BWA_HOME)

$(SNAP):
	$(MAKE) -C $(BWA_HOME)

include ../mason.mk
include ../bowtie2ts.mk
include ../bwamemts.mk
include ../snapts.mk

.PHONY: clean
clean:
	rm -f *.fq.gz
	rm -rf *.ts *.pred

.PHONY: out_bt2 out_snap out_bwamem

out_bt2: \
    r0_bt2s_mason_ill_100_50M.out \
    r0_bt2s_mason_ill_250_50M.out \
    r12_bt2s100_mason_ill_100_50M.out \
    r12_bt2s250_mason_ill_250_50M.out

out_snap: \
    r0_snap_mason_ill_100_50M.out \
    r0_snap_mason_ill_250_50M.out \
    r12_snap100_mason_ill_100_50M.out \
    r12_snap250_mason_ill_250_50M.out

out_bwamem: \
    r0_bwamem_mason_ill_100_50M.out \
    r0_bwamem_mason_ill_250_50M.out \
    r12_bwamem100_mason_ill_100_50M.out \
    r12_bwamem250_mason_ill_250_50M.out

.PHONY: reads
reads: \
    r0_mason_ill_100_50M.fq.gz \
    r0_mason_ill_250_50M.fq.gz \
    r1_mason_ill_100_50M.fq.gz \
    r1_mason_ill_250_50M.fq.gz

$(eval $(call bt2ts,bt2s,1,1,,,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2s100,1,1,,-I 200 -X 400,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2s250,1,1,,-I 500 -X 1000,$(HG_FASTA),$(HG_INDEX)))

$(eval $(call bwamemts,bwamem,1,1,,,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bwamemts,bwamem100,1,1,,-I 300$(,)25,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bwamemts,bwamem250,1,1,,-I 750$(,)62,$(HG_FASTA),$(HG_INDEX)))

$(eval $(call snapts,snap,1,1,,,,,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call snapts,snap100,1,1,,,,-s 100 300,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call snapts,snap250,1,1,,,,-s 250 750,$(HG_FASTA),$(HG_INDEX)))

$(eval $(call mason_ill_unp_reads,ill_100_50M,$(HG_FASTA),100,50000000,7724))
$(eval $(call mason_ill_unp_reads,ill_250_50M,$(HG_FASTA),250,50000000,7724))

$(eval $(call mason_ill_pair_reads,ill_100_50M,$(HG_FASTA),100,50000000,300,100,7724))
$(eval $(call mason_ill_pair_reads,ill_250_50M,$(HG_FASTA),250,50000000,750,250,7724))
