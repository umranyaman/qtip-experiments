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

.PHONY: outs
outs: outs_bt2

.PHONY: clean
clean:
	rm -f *.fq.gz
	rm -rf *.out

.PHONY: outs
outs: outs_unpaired outs_paired

.PHONY: outs_unpaired
outs_unpaired: \
    r0_bt2s_mason_ill_50.out/DONE \
    r0_bt2s_mason_ill_100.out/DONE \
    r0_bt2s_mason_ill_150.out/DONE \
    r0_bt2s_mason_ill_250.out/DONE \
    r0_bt2s_mason_ill_500.out/DONE \
    r0_bt2sl_mason_ill_50.out/DONE \
    r0_bt2sl_mason_ill_100.out/DONE \
    r0_bt2sl_mason_ill_150.out/DONE \
    r0_bt2sl_mason_ill_250.out/DONE \
    r0_bt2sl_mason_ill_500.out/DONE

.PHONY: outs_paired
outs_paired: \
    r12_bt2s50_mason_ill_50.out/DONE \
    r12_bt2s100_mason_ill_100.out/DONE \
    r12_bt2s150_mason_ill_150.out/DONE \
    r12_bt2s250_mason_ill_250.out/DONE \
    r12_bt2s500_mason_ill_500.out/DONE \
    r12_bt2sl50_mason_ill_50.out/DONE \
    r12_bt2sl100_mason_ill_100.out/DONE \
    r12_bt2sl150_mason_ill_150.out/DONE \
    r12_bt2sl250_mason_ill_250.out/DONE \
    r12_bt2sl500_mason_ill_500.out/DONE

.PHONY: reads
reads: \
    r0_mason_ill_50.fq.gz  \
    r0_mason_ill_100.fq.gz \
    r0_mason_ill_150.fq.gz \
    r0_mason_ill_250.fq.gz \
    r0_mason_ill_500.fq.gz \
    r1_mason_ill_50.fq.gz  \
    r1_mason_ill_100.fq.gz \
    r1_mason_ill_150.fq.gz \
    r1_mason_ill_250.fq.gz \
    r1_mason_ill_500.fq.gz

SUBSAMPLING_ARGS=--subsampling-series 0.015,0.03 --trials=5 --predict-for-training --keep-intermediates --assess-accuracy

#
# Make bowtie2 rule.  Use --sensitive, hg19, and prefix rules with "bt2s"
#
$(eval $(call bt2ts,bt2s,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),--sensitive,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2sl,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),--sensitive --local,$(HG_FASTA),$(HG_INDEX),out))

$(eval $(call bt2ts,bt2s50,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 100 -X 200 --sensitive,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2sl50,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 100 -X 200 --sensitive --local,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2s100,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 200 -X 400 --sensitive,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2sl100,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 200 -X 400 --sensitive --local,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2s150,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 300 -X 600 --sensitive,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2sl150,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 300 -X 600 --sensitive --local,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2s250,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 500 -X 1000 --sensitive,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2sl250,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 500 -X 1000 --sensitive --local,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2s500,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 1000 -X 2000 --sensitive,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2sl500,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 1000 -X 2000 --sensitive --local,$(HG_FASTA),$(HG_INDEX),out))

$(eval $(call mason_ill_unp_reads,ill_50,$(HG_FASTA),50,$(NREADS),7724))
$(eval $(call mason_ill_unp_reads,ill_100,$(HG_FASTA),100,$(NREADS),7724))
$(eval $(call mason_ill_unp_reads,ill_150,$(HG_FASTA),150,$(NREADS),7724))
$(eval $(call mason_ill_unp_reads,ill_250,$(HG_FASTA),250,$(NREADS),7724))
$(eval $(call mason_ill_unp_reads,ill_500,$(HG_FASTA),500,$(NREADS),7724))

$(eval $(call mason_ill_pair_reads,ill_50,$(HG_FASTA),50,$(NREADS),150,50,7724))
$(eval $(call mason_ill_pair_reads,ill_100,$(HG_FASTA),100,$(NREADS),300,100,7724))
$(eval $(call mason_ill_pair_reads,ill_150,$(HG_FASTA),150,$(NREADS),450,150,7724))
$(eval $(call mason_ill_pair_reads,ill_250,$(HG_FASTA),250,$(NREADS),750,250,7724))
$(eval $(call mason_ill_pair_reads,ill_500,$(HG_FASTA),500,$(NREADS),1500,500,7724))
