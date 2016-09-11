#!gmake

# Generate simulated FASTQ files for Illumina-like reads, both
# unpaired and paired-end.  We generate several files, each for a
# different read length from 50 up to 500 nt.  Uses the Mason
# simulator.

# Prerequisites:
# 1. QTIP_EXPERIMENTS_HOME (root of repo) must be set
# 2. $QTIP_EXPERIMENTS_HOME/experiments/refs must have relevant references and indexes
# 3. $QTIP_EXPERIMENTS_HOME/software/mason/mason must have been built

HG_FASTA_NAME=hg19.fa
MM_FASTA_NAME=mm10.fa
ZM_FASTA_NAME=zm_AGPv3.fa
HG_INDEX_BASE=$(HG_FASTA_NAME)
MM_INDEX_BASE=$(MM_FASTA_NAME)
ZM_INDEX_BASE=$(ZM_FASTA_NAME)

ifeq (,$(QTIP_EXPERIMENTS_HOME))
$(error "Set QTIP_EXPERIMENTS_HOME, where qtip-experiments working copy is located")
endif
TS_REFS=$(QTIP_EXPERIMENTS_HOME)/experiments/refs
TS_INDEXES=$(QTIP_EXPERIMENTS_HOME)/experiments/refs

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

BOWTIE2_HOME=$(QTIP_EXPERIMENTS_HOME)/software/bowtie2
BOWTIE2=$(BOWTIE2_HOME)/bowtie2
BT2_ARGS=

BWA_HOME=$(QTIP_EXPERIMENTS_HOME)/software/bwa
BWA=$(BWA_HOME)/bwa
BWA_ARGS=

SNAP_HOME=$(QTIP_EXPERIMENTS_HOME)/software/snap
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

NREADS=4000000
NCORES=1

.PHONY: clean
clean:
	rm -f *.fq.gz
	rm -rf *.out

.PHONY: outs
outs: outs_bt2 outs_bwamem outs_snap

.PHONY: outs_bt2
outs_bt2: \
    r0_bt2s_hg_mason_ill_hg_100.out/DONE \
    r0_bt2s_hg_mason_ill_hg_250.out/DONE \
    r0_bt2s_mm_mason_ill_mm_100.out/DONE \
    r0_bt2s_mm_mason_ill_mm_250.out/DONE \
    r0_bt2s_zm_mason_ill_zm_100.out/DONE \
    r0_bt2s_zm_mason_ill_zm_250.out/DONE \
    r12_bt2s100_hg_mason_ill_hg_100.out/DONE \
    r12_bt2s250_hg_mason_ill_hg_250.out/DONE \
    r12_bt2s100_mm_mason_ill_mm_100.out/DONE \
    r12_bt2s250_mm_mason_ill_mm_250.out/DONE \
    r12_bt2s100_zm_mason_ill_zm_100.out/DONE \
    r12_bt2s250_zm_mason_ill_zm_250.out/DONE

.PHONY: outs_bwamem
outs_bwamem: \
    r0_bwamem_hg_mason_ill_hg_100.out/DONE \
    r0_bwamem_hg_mason_ill_hg_250.out/DONE \
    r0_bwamem_mm_mason_ill_mm_100.out/DONE \
    r0_bwamem_mm_mason_ill_mm_250.out/DONE \
    r0_bwamem_zm_mason_ill_zm_100.out/DONE \
    r0_bwamem_zm_mason_ill_zm_250.out/DONE \
    r12_bwamem100_hg_mason_ill_hg_100.out/DONE \
    r12_bwamem250_hg_mason_ill_hg_250.out/DONE \
    r12_bwamem100_mm_mason_ill_mm_100.out/DONE \
    r12_bwamem250_mm_mason_ill_mm_250.out/DONE \
    r12_bwamem100_zm_mason_ill_zm_100.out/DONE \
    r12_bwamem250_zm_mason_ill_zm_250.out/DONE

.PHONY: outs_snap
outs_snap: \
    r0_snap_hg_mason_ill_hg_100.out/DONE \
    r0_snap_hg_mason_ill_hg_250.out/DONE \
    r0_snap_mm_mason_ill_mm_100.out/DONE \
    r0_snap_mm_mason_ill_mm_250.out/DONE \
    r0_snap_zm_mason_ill_zm_100.out/DONE \
    r0_snap_zm_mason_ill_zm_250.out/DONE \
    r12_snap100_hg_mason_ill_hg_100.out/DONE \
    r12_snap250_hg_mason_ill_hg_250.out/DONE \
    r12_snap100_mm_mason_ill_mm_100.out/DONE \
    r12_snap250_mm_mason_ill_mm_250.out/DONE \
    r12_snap100_zm_mason_ill_zm_100.out/DONE \
    r12_snap250_zm_mason_ill_zm_250.out/DONE

.PHONY: reads
reads: \
    r0_mason_ill_hg_100.fq.gz \
    r0_mason_ill_hg_250.fq.gz \
    r0_mason_ill_mm_100.fq.gz \
    r0_mason_ill_mm_250.fq.gz \
    r0_mason_ill_zm_100.fq.gz \
    r0_mason_ill_zm_250.fq.gz \
    r1_mason_ill_hg_100.fq.gz \
    r1_mason_ill_hg_250.fq.gz \
    r1_mason_ill_mm_100.fq.gz \
    r1_mason_ill_mm_250.fq.gz \
    r1_mason_ill_zm_100.fq.gz \
    r1_mason_ill_zm_250.fq.gz

MK_QTIP_ARGS=--trials=10 --predict-for-training --keep-intermediates --assess-accuracy

, := ,

$(eval $(call bt2ts,bt2s_hg,$(MK_QTIP_ARGS),--sensitive,$(FA_HG),$(HG_INDEX),out,$(NCORES)))
$(eval $(call bt2ts,bt2s_mm,$(MK_QTIP_ARGS),--sensitive,$(FA_MM),$(MM_INDEX),out,$(NCORES)))
$(eval $(call bt2ts,bt2s_zm,$(MK_QTIP_ARGS),--sensitive,$(FA_ZM),$(ZM_INDEX),out,$(NCORES)))

$(eval $(call bwamemts,bwamem_hg,$(MK_QTIP_ARGS),,$(FA_HG),$(HG_INDEX),out,$(NCORES)))
$(eval $(call bwamemts,bwamem_mm,$(MK_QTIP_ARGS),,$(FA_MM),$(MM_INDEX),out,$(NCORES)))
$(eval $(call bwamemts,bwamem_zm,$(MK_QTIP_ARGS),,$(FA_ZM),$(ZM_INDEX),out,$(NCORES)))

$(eval $(call snapts,snap_hg,$(MK_QTIP_ARGS),,,,$(FA_HG),$(HG_INDEX),out,$(NCORES)))
$(eval $(call snapts,snap_mm,$(MK_QTIP_ARGS),,,,$(FA_MM),$(MM_INDEX),out,$(NCORES)))
$(eval $(call snapts,snap_zm,$(MK_QTIP_ARGS),,,,$(FA_ZM),$(ZM_INDEX),out,$(NCORES)))

$(eval $(call bt2ts,bt2s100_hg,$(MK_QTIP_ARGS),-I 200 -X 400 --sensitive,$(FA_HG),$(HG_INDEX),out,$(NCORES)))
$(eval $(call bt2ts,bt2s100_mm,$(MK_QTIP_ARGS),-I 200 -X 400 --sensitive,$(FA_MM),$(MM_INDEX),out,$(NCORES)))
$(eval $(call bt2ts,bt2s100_zm,$(MK_QTIP_ARGS),-I 200 -X 400 --sensitive,$(FA_ZM),$(ZM_INDEX),out,$(NCORES)))

$(eval $(call bwamemts,bwamem100_hg,$(MK_QTIP_ARGS),-I 300$(,)25,$(FA_HG),$(HG_INDEX),out,$(NCORES)))
$(eval $(call bwamemts,bwamem100_mm,$(MK_QTIP_ARGS),-I 300$(,)25,$(FA_MM),$(MM_INDEX),out,$(NCORES)))
$(eval $(call bwamemts,bwamem100_zm,$(MK_QTIP_ARGS),-I 300$(,)25,$(FA_ZM),$(ZM_INDEX),out,$(NCORES)))

$(eval $(call snapts,snap100_hg,$(MK_QTIP_ARGS),,,-s 100 300,$(FA_HG),$(HG_INDEX),out,$(NCORES)))
$(eval $(call snapts,snap100_mm,$(MK_QTIP_ARGS),,,-s 100 300,$(FA_MM),$(MM_INDEX),out,$(NCORES)))
$(eval $(call snapts,snap100_zm,$(MK_QTIP_ARGS),,,-s 100 300,$(FA_ZM),$(ZM_INDEX),out,$(NCORES)))

$(eval $(call bt2ts,bt2s250_hg,$(MK_QTIP_ARGS),-I 500 -X 1000 --sensitive,$(FA_HG),$(HG_INDEX),out,$(NCORES)))
$(eval $(call bt2ts,bt2s250_mm,$(MK_QTIP_ARGS),-I 500 -X 1000 --sensitive,$(FA_MM),$(MM_INDEX),out,$(NCORES)))
$(eval $(call bt2ts,bt2s250_zm,$(MK_QTIP_ARGS),-I 500 -X 1000 --sensitive,$(FA_ZM),$(ZM_INDEX),out,$(NCORES)))

$(eval $(call bwamemts,bwamem250_hg,$(MK_QTIP_ARGS),-I 750$(,)62,$(FA_HG),$(HG_INDEX),out,$(NCORES)))
$(eval $(call bwamemts,bwamem250_mm,$(MK_QTIP_ARGS),-I 750$(,)62,$(FA_MM),$(MM_INDEX),out,$(NCORES)))
$(eval $(call bwamemts,bwamem250_zm,$(MK_QTIP_ARGS),-I 750$(,)62,$(FA_ZM),$(ZM_INDEX),out,$(NCORES)))

$(eval $(call snapts,snap250_hg,$(MK_QTIP_ARGS),,,-s 250 750,$(FA_HG),$(HG_INDEX),out,$(NCORES)))
$(eval $(call snapts,snap250_mm,$(MK_QTIP_ARGS),,,-s 250 750,$(FA_MM),$(MM_INDEX),out,$(NCORES)))
$(eval $(call snapts,snap250_zm,$(MK_QTIP_ARGS),,,-s 250 750,$(FA_ZM),$(ZM_INDEX),out,$(NCORES)))

$(eval $(call mason_ill_unp_reads,ill_hg_100,$(FA_HG),100,$(NREADS),7724))
$(eval $(call mason_ill_unp_reads,ill_hg_250,$(FA_HG),250,$(NREADS),7724))

$(eval $(call mason_ill_unp_reads,ill_mm_100,$(FA_MM),100,$(NREADS),7724))
$(eval $(call mason_ill_unp_reads,ill_mm_250,$(FA_MM),250,$(NREADS),7724))

$(eval $(call mason_ill_unp_reads,ill_zm_100,$(FA_ZM),100,$(NREADS),7724))
$(eval $(call mason_ill_unp_reads,ill_zm_250,$(FA_ZM),250,$(NREADS),7724))

$(eval $(call mason_ill_pair_reads,ill_hg_100,$(FA_HG),100,$(NREADS),300,100,7724))
$(eval $(call mason_ill_pair_reads,ill_hg_250,$(FA_HG),250,$(NREADS),750,250,7724))

$(eval $(call mason_ill_pair_reads,ill_mm_100,$(FA_MM),100,$(NREADS),300,100,7724))
$(eval $(call mason_ill_pair_reads,ill_mm_250,$(FA_MM),250,$(NREADS),750,250,7724))

$(eval $(call mason_ill_pair_reads,ill_zm_100,$(FA_ZM),100,$(NREADS),300,100,7724))
$(eval $(call mason_ill_pair_reads,ill_zm_250,$(FA_ZM),250,$(NREADS),750,250,7724))
