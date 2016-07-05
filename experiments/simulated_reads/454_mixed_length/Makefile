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

.PHONY: all
all: outs

BOWTIE2_HOME=$(QSIM_EXPERIMENTS_HOME)/software/bowtie2
BOWTIE2=$(BOWTIE2_HOME)/bowtie2
BT2_ARGS=

$(BOWTIE2):
	$(MAKE) -C $(BOWTIE2_HOME)

include ../mason.mk
include ../bowtie2ts.mk

NREADS=4000000
NCORES=8

.PHONY: clean
clean:
	rm -f *.fq.gz
	rm -rf *.out

.PHONY: outs
outs: outs_bt2

.PHONY: outs_bt2
outs_bt2: \
    r0_bt2s_mason_fff_50to500.out/DONE \
    r0_bt2sl_mason_fff_50to500.out/DONE

.PHONY: reads
reads: \
    r0_mason_fff_50to500.fq.gz

SUBSAMPLING_ARGS=--subsampling-series 0.015,0.03 --trials=5 --predict-for-training --keep-intermediates --assess-accuracy

#
# Make bowtie2 rule.  Use --sensitive, hg19, and prefix rules with "bt2s"
#
$(eval $(call bt2ts,bt2s,$(SUBSAMPLING_ARGS),--sensitive,$(HG_FASTA),$(HG_INDEX),out,$(NCORES)))
$(eval $(call bt2ts,bt2sl,$(SUBSAMPLING_ARGS),--sensitive --local,$(HG_FASTA),$(HG_INDEX),out,$(NCORES)))

$(eval $(call mason_fff_unp_reads,fff_50to500,$(HG_FASTA),275,225,$(NREADS),7724))
