#!gmake

# Generate simulated FASTQ files for Illumina-like reads, both
# unpaired and paired-end.  We generate several files, each for a
# different read length from 50 up to 500 nt.  Uses the Mason
# simulator.

# Prerequisites:
# 1. QSIM_EXPERIMENTS_HOME (root of repo) must be set
# 2. $QSIM_EXPERIMENTS_HOME/experiments/refs must have relevant references and indexes
# 3. $QSIM_EXPERIMENTS_HOME/software/mason/mason must have been built

HG_FASTA_NAME=hg19.fa
HG_INDEX_BASE=$(HG_FASTA_NAME)

ifeq (,$(QSIM_EXPERIMENTS_HOME))
$(error "Set QSIM_EXPERIMENTS_HOME, where qsim-experiments working copy is located")
endif
TS_REFS=$(QSIM_EXPERIMENTS_HOME)/experiments/refs
TS_INDEXES=$(QSIM_EXPERIMENTS_HOME)/experiments/refs

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
all: outs

BOWTIE2_HOME=$(QSIM_EXPERIMENTS_HOME)/software/bowtie2
BOWTIE2=$(BOWTIE2_HOME)/bowtie2
BT2_ARGS=

BWA_HOME=$(QSIM_EXPERIMENTS_HOME)/software/bwa
BWA=$(BWA_HOME)/bwa
BWA_ARGS=

SNAP_HOME=$(QSIM_EXPERIMENTS_HOME)/software/snap
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

.PHONY: clean
clean:
	rm -f *.fq.gz
	rm -rf *.out

.PHONY: outs
outs: outs_bt2 outs_bwamem outs_snap

.PHONY: outs_bt2
outs_bt2: \
    r0_bt2s_mason_ill_100.out/DONE \
    r0_bt2s_mason_ill_250.out/DONE \
    r12_bt2s100_mason_ill_100.out/DONE \
    r12_bt2s250_mason_ill_250.out/DONE \
    r0_bt2sl_mason_ill_100.out/DONE \
    r0_bt2sl_mason_ill_250.out/DONE \
    r12_bt2sl100_mason_ill_100.out/DONE \
    r12_bt2sl250_mason_ill_250.out/DONE

.PHONY: outs_bwamem
outs_bwamem: \
    r0_bwamem_mason_ill_100.out/DONE \
    r0_bwamem_mason_ill_250.out/DONE \
    r12_bwamem100_mason_ill_100.out/DONE \
    r12_bwamem250_mason_ill_250.out/DONE

.PHONY: outs_snap
outs_snap: \
    r0_snap_mason_ill_100.out/DONE \
    r0_snap_mason_ill_250.out/DONE \
    r12_snap100_mason_ill_100.out/DONE \
    r12_snap250_mason_ill_250.out/DONE

.PHONY: reads
reads: \
    r0_mason_ill_100.fq.gz \
    r0_mason_ill_250.fq.gz \
    r1_mason_ill_100.fq.gz \
    r1_mason_ill_250.fq.gz

SUBSAMPLING_ARGS=--subsampling-series 0.015,0.03 --trials=5 --predict-for-training --keep-intermediates --assess-accuracy

#
# Make bowtie2 rule.  Use --sensitive, hg19, and prefix rules with "bt2s"
#
$(eval $(call bt2ts,bt2s,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),--sensitive,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2sl,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),--sensitive --local,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2s100,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 200 -X 400 --sensitive,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2sl100,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 200 -X 400 --sensitive --local,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2s250,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 500 -X 1000 --sensitive,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2sl250,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 500 -X 1000 --sensitive --local,$(HG_FASTA),$(HG_INDEX),out))

, := ,

$(eval $(call bwamemts,bwamem,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bwamemts,bwamem100,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 300$(,)25,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bwamemts,bwamem250,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 750$(,)62,$(HG_FASTA),$(HG_INDEX),out))

$(eval $(call snapts,snap,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),,,,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call snapts,snap100,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),,,-s 100 300,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call snapts,snap250,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),,,-s 250 750,$(HG_FASTA),$(HG_INDEX),out))

$(eval $(call mason_ill_unp_reads,ill_100,$(HG_FASTA),100,$(NREADS),7724))
$(eval $(call mason_ill_unp_reads,ill_250,$(HG_FASTA),250,$(NREADS),7724))

$(eval $(call mason_ill_pair_reads,ill_100,$(HG_FASTA),100,$(NREADS),300,100,7724))
$(eval $(call mason_ill_pair_reads,ill_250,$(HG_FASTA),250,$(NREADS),750,250,7724))
