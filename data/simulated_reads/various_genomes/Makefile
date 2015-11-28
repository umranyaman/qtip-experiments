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
MM_FASTA_NAME=mm10.fa
ZM_FASTA_NAME=zm_AGPv3.fa
HG_INDEX_BASE=$(HG_FASTA_NAME)
MM_INDEX_BASE=$(MM_FASTA_NAME)
ZM_INDEX_BASE=$(ZM_FASTA_NAME)

ifeq (,$(TS_HOME))
$(error "Set TS_HOME, where ts repo root is located")
endif
ifeq (,$(TS_REFS))
$(error "Set TS_REFS, where $(REF_FASTA) is located")
endif
ifeq (,$(TS_INDEXES))
$(error "Set TS_INDEXES, where ts repo root is located")
endif

FA_HG=$(TS_REFS)/$(HG_FASTA_NAME)
FA_MM=$(TS_REFS)/$(MM_FASTA_NAME)
FA_ZM=$(TS_REFS)/$(ZM_FASTA_NAME)

HG_INDEX=$(TS_INDEXES)/$(HG_INDEX_BASE)
MM_INDEX=$(TS_INDEXES)/$(MM_INDEX_BASE)
ZM_INDEX=$(TS_INDEXES)/$(ZM_INDEX_BASE)

ifeq ($(wildcard $(FA_HG)),) 
$(error "$(FA_HG) doesn't exist")
endif
#ifeq ($(wildcard $(FA_MM)),) 
#$(error "$(FA_MM) doesn't exist")
#endif
ifeq ($(wildcard $(FA_ZM)),) 
$(error "$(FA_ZM) doesn't exist")
endif
ifeq ($(wildcard $(HG_INDEX).1.bt2),) 
$(error "$(HG_INDEX).1.bt2 doesn't exist")
endif
#ifeq ($(wildcard $(MM_INDEX).1.bt2),) 
#$(error "$(MM_INDEX).1.bt2 doesn't exist")
#endif
ifeq ($(wildcard $(ZM_INDEX).1.bt2),) 
$(error "$(ZM_INDEX).1.bt2 doesn't exist")
endif

.PHONY: all
all: outs

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
	$(MAKE) -C $(SNAP_HOME)


include ../mason.mk
include ../bowtie2ts.mk
include ../bwamemts.mk
include ../snapts.mk

.PHONY: clean
clean:
	rm -f *.fq.gz
	rm -f *.mkb
	rm -rf *.out

.PHONY: outs
outs: outs_bt2 outs_bwamem

.PHONY: outs_bt2
outs_bt2: \
    r0_bt2s_hg_mason_ill_hg_100_100k.out \
    r0_bt2s_hg_mason_ill_hg_250_100k.out \
    r0_bt2s_mm_mason_ill_mm_100_100k.out \
    r0_bt2s_mm_mason_ill_mm_250_100k.out \
    r0_bt2s_zm_mason_ill_zm_100_100k.out \
    r0_bt2s_zm_mason_ill_zm_250_100k.out \
    r12_bt2s100_hg_mason_ill_hg_100_100k.out \
    r12_bt2s250_hg_mason_ill_hg_250_100k.out \
    r12_bt2s100_mm_mason_ill_mm_100_100k.out \
    r12_bt2s250_mm_mason_ill_mm_250_100k.out \
    r12_bt2s100_zm_mason_ill_zm_100_100k.out \
    r12_bt2s250_zm_mason_ill_zm_250_100k.out

.PHONY: outs_bwamem
outs_bwamem: \
    r0_bwamem_hg_mason_ill_hg_100_100k.out \
    r0_bwamem_hg_mason_ill_hg_250_100k.out \
    r0_bwamem_mm_mason_ill_mm_100_100k.out \
    r0_bwamem_mm_mason_ill_mm_250_100k.out \
    r0_bwamem_zm_mason_ill_zm_100_100k.out \
    r0_bwamem_zm_mason_ill_zm_250_100k.out \
    r12_bwamem100_hg_mason_ill_hg_100_100k.out \
    r12_bwamem250_hg_mason_ill_hg_250_100k.out \
    r12_bwamem100_mm_mason_ill_mm_100_100k.out \
    r12_bwamem250_mm_mason_ill_mm_250_100k.out \
    r12_bwamem100_zm_mason_ill_zm_100_100k.out \
    r12_bwamem250_zm_mason_ill_zm_250_100k.out

SUBSAMPLING_ARGS=--subsampling-series 0.05,0.1,0.2,0.3,0.4,0.5,1.0 --trials=5 --predict-for-training --keep-intermediates --assess-accuracy

, := ,

$(eval $(call bt2ts,bt2s_hg,100000,100000,$(SUBSAMPLING_ARGS),--sensitive,$(FA_HG),$(HG_INDEX)))
$(eval $(call bt2ts,bt2s_mm,100000,100000,$(SUBSAMPLING_ARGS),--sensitive,$(FA_MM),$(MM_INDEX)))
$(eval $(call bt2ts,bt2s_zm,100000,100000,$(SUBSAMPLING_ARGS),--sensitive,$(FA_ZM),$(ZM_INDEX)))

$(eval $(call bwamemts,bwamem_hg,100000,100000,$(SUBSAMPLING_ARGS),,$(FA_HG),$(HG_INDEX)))
$(eval $(call bwamemts,bwamem_mm,100000,100000,$(SUBSAMPLING_ARGS),,$(FA_MM),$(MM_INDEX)))
$(eval $(call bwamemts,bwamem_zm,100000,100000,$(SUBSAMPLING_ARGS),,$(FA_ZM),$(ZM_INDEX)))

$(eval $(call bt2ts,bt2s100_hg,100000,100000,$(SUBSAMPLING_ARGS),-I 200 -X 400 --sensitive,$(FA_HG),$(HG_INDEX)))
$(eval $(call bt2ts,bt2s100_mm,100000,100000,$(SUBSAMPLING_ARGS),-I 200 -X 400 --sensitive,$(FA_MM),$(MM_INDEX)))
$(eval $(call bt2ts,bt2s100_zm,100000,100000,$(SUBSAMPLING_ARGS),-I 200 -X 400 --sensitive,$(FA_ZM),$(ZM_INDEX)))

$(eval $(call bwamemts,bwamem100_hg,100000,100000,$(SUBSAMPLING_ARGS),-I 300$(,)25,$(FA_HG),$(HG_INDEX)))
$(eval $(call bwamemts,bwamem100_mm,100000,100000,$(SUBSAMPLING_ARGS),-I 300$(,)25,$(FA_MM),$(MM_INDEX)))
$(eval $(call bwamemts,bwamem100_zm,100000,100000,$(SUBSAMPLING_ARGS),-I 300$(,)25,$(FA_ZM),$(ZM_INDEX)))

$(eval $(call bt2ts,bt2s250_hg,100000,100000,$(SUBSAMPLING_ARGS),-I 500 -X 1000 --sensitive,$(FA_HG),$(HG_INDEX)))
$(eval $(call bt2ts,bt2s250_mm,100000,100000,$(SUBSAMPLING_ARGS),-I 500 -X 1000 --sensitive,$(FA_MM),$(MM_INDEX)))
$(eval $(call bt2ts,bt2s250_zm,100000,100000,$(SUBSAMPLING_ARGS),-I 500 -X 1000 --sensitive,$(FA_ZM),$(ZM_INDEX)))

$(eval $(call bwamemts,bwamem250_hg,100000,100000,$(SUBSAMPLING_ARGS),-I 750$(,)62,$(FA_HG),$(HG_INDEX)))
$(eval $(call bwamemts,bwamem250_mm,100000,100000,$(SUBSAMPLING_ARGS),-I 750$(,)62,$(FA_MM),$(MM_INDEX)))
$(eval $(call bwamemts,bwamem250_zm,100000,100000,$(SUBSAMPLING_ARGS),-I 750$(,)62,$(FA_ZM),$(ZM_INDEX)))

$(eval $(call mason_ill_unp_reads,ill_hg_100_100k,$(FA_HG),100,100000,7724))
$(eval $(call mason_ill_unp_reads,ill_hg_250_100k,$(FA_HG),250,100000,7724))

$(eval $(call mason_ill_unp_reads,ill_mm_100_100k,$(FA_MM),100,100000,7724))
$(eval $(call mason_ill_unp_reads,ill_mm_250_100k,$(FA_MM),250,100000,7724))

$(eval $(call mason_ill_unp_reads,ill_zm_100_100k,$(FA_ZM),100,100000,7724))
$(eval $(call mason_ill_unp_reads,ill_zm_250_100k,$(FA_ZM),250,100000,7724))

$(eval $(call mason_ill_pair_reads,ill_hg_100_100k,$(FA_HG),100,100000,300,100,7724))
$(eval $(call mason_ill_pair_reads,ill_hg_250_100k,$(FA_HG),250,100000,750,250,7724))

$(eval $(call mason_ill_pair_reads,ill_mm_100_100k,$(FA_MM),100,100000,300,100,7724))
$(eval $(call mason_ill_pair_reads,ill_mm_250_100k,$(FA_MM),250,100000,750,250,7724))

$(eval $(call mason_ill_pair_reads,ill_zm_100_100k,$(FA_ZM),100,100000,300,100,7724))
$(eval $(call mason_ill_pair_reads,ill_zm_250_100k,$(FA_ZM),250,100000,750,250,7724))
