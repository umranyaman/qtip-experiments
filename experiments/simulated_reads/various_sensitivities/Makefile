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
NREADS_NM=4M

.PHONY: clean
clean:
	rm -f *.fq.gz
	rm -rf *.out

.PHONY: outs
outs: outs_unpaired outs_paired

.PHONY: outs_unpaired
outs_unpaired: \
    r0_bt2vf_mason_ill_100_$(NREADS_NM).out/DONE \
    r0_bt2vf_mason_ill_250_$(NREADS_NM).out/DONE \
    r0_bt2f_mason_ill_100_$(NREADS_NM).out/DONE \
    r0_bt2f_mason_ill_250_$(NREADS_NM).out/DONE \
    r0_bt2s_mason_ill_100_$(NREADS_NM).out/DONE \
    r0_bt2s_mason_ill_250_$(NREADS_NM).out/DONE \
    r0_bt2vs_mason_ill_100_$(NREADS_NM).out/DONE \
    r0_bt2vs_mason_ill_250_$(NREADS_NM).out/DONE \
    r0_bt2vfl_mason_ill_100_$(NREADS_NM).out/DONE \
    r0_bt2vfl_mason_ill_250_$(NREADS_NM).out/DONE \
    r0_bt2fl_mason_ill_100_$(NREADS_NM).out/DONE \
    r0_bt2fl_mason_ill_250_$(NREADS_NM).out/DONE \
    r0_bt2sl_mason_ill_100_$(NREADS_NM).out/DONE \
    r0_bt2sl_mason_ill_250_$(NREADS_NM).out/DONE \
    r0_bt2vsl_mason_ill_100_$(NREADS_NM).out/DONE \
    r0_bt2vsl_mason_ill_250_$(NREADS_NM).out/DONE

.PHONY: outs_paired
outs_paired: \
    r12_bt2vf100_mason_ill_100_$(NREADS_NM).out/DONE \
    r12_bt2vf250_mason_ill_250_$(NREADS_NM).out/DONE \
    r12_bt2f100_mason_ill_100_$(NREADS_NM).out/DONE \
    r12_bt2f250_mason_ill_250_$(NREADS_NM).out/DONE \
    r12_bt2s100_mason_ill_100_$(NREADS_NM).out/DONE \
    r12_bt2s250_mason_ill_250_$(NREADS_NM).out/DONE \
    r12_bt2vs100_mason_ill_100_$(NREADS_NM).out/DONE \
    r12_bt2vs250_mason_ill_250_$(NREADS_NM).out/DONE \
    r12_bt2vfl100_mason_ill_100_$(NREADS_NM).out/DONE \
    r12_bt2vfl250_mason_ill_250_$(NREADS_NM).out/DONE \
    r12_bt2fl100_mason_ill_100_$(NREADS_NM).out/DONE \
    r12_bt2fl250_mason_ill_250_$(NREADS_NM).out/DONE \
    r12_bt2sl100_mason_ill_100_$(NREADS_NM).out/DONE \
    r12_bt2sl250_mason_ill_250_$(NREADS_NM).out/DONE \
    r12_bt2vsl100_mason_ill_100_$(NREADS_NM).out/DONE \
    r12_bt2vsl250_mason_ill_250_$(NREADS_NM).out/DONE

.PHONY: reads
reads: \
    r0_mason_ill_100_$(NREADS_NM).fq.gz \
    r0_mason_ill_250_$(NREADS_NM).fq.gz \
    r1_mason_ill_100_$(NREADS_NM).fq.gz \
    r1_mason_ill_250_$(NREADS_NM).fq.gz

SUBSAMPLING_ARGS=--subsampling-series 0.015,0.03 --trials=5 --predict-for-training --keep-intermediates --assess-accuracy

#
# Make bowtie2 rule.  Use --sensitive, hg19, and prefix rules with "bt2s"
#
$(eval $(call bt2ts,bt2vf,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),--very-fast,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2f,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),--fast,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2s,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),--sensitive,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2vs,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),--very-sensitive,$(HG_FASTA),$(HG_INDEX),out))

$(eval $(call bt2ts,bt2vfl,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),--very-fast-local,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2fl,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),--fast-local,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2sl,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),--sensitive-local,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2vsl,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),--very-sensitive-local,$(HG_FASTA),$(HG_INDEX),out))

$(eval $(call bt2ts,bt2vf100,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 200 -X 400 --very-fast,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2vf250,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 500 -X 1000 --very-fast,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2f100,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 200 -X 400 --fast,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2f250,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 500 -X 1000 --fast,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2s100,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 200 -X 400 --sensitive,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2s250,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 500 -X 1000 --sensitive,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2vs100,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 200 -X 400 --very-sensitive,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2vs250,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 500 -X 1000 --very-sensitive,$(HG_FASTA),$(HG_INDEX),out))

$(eval $(call bt2ts,bt2vfl100,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 200 -X 400 --very-fast-local,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2vfl250,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 500 -X 1000 --very-fast-local,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2fl100,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 200 -X 400 --fast-local,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2fl250,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 500 -X 1000 --fast-local,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2sl100,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 200 -X 400 --sensitive-local,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2sl250,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 500 -X 1000 --sensitive-local,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2vsl100,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 200 -X 400 --very-sensitive-local,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2vsl250,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 500 -X 1000 --very-sensitive-local,$(HG_FASTA),$(HG_INDEX),out))

$(eval $(call mason_ill_unp_reads,ill_100_$(NREADS_NM),$(HG_FASTA),100,$(NREADS),7724))
$(eval $(call mason_ill_unp_reads,ill_250_$(NREADS_NM),$(HG_FASTA),250,$(NREADS),7724))

$(eval $(call mason_ill_pair_reads,ill_100_$(NREADS_NM),$(HG_FASTA),100,$(NREADS),300,100,7724))
$(eval $(call mason_ill_pair_reads,ill_250_$(NREADS_NM),$(HG_FASTA),250,$(NREADS),750,250,7724))
