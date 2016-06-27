#!gmake

# Generate simulated FASTQ files for Illumina-like reads, both
# unpaired and paired-end.  We generate several files, each for a
# different read length from 50 up to 500 nt.  Uses the Mason
# simulator.

# Prerequisites:
# 1. QSIM_EXPERIMENTS_HOME (root of repo) must be set
# 2. $QSIM_EXPERIMENTS_HOME/experiments/refs must have relevant references and indexes
# 3. $QSIM_EXPERIMENTS_HOME/software/wgsim/wgsim must have been built
# 4. $QSIM_EXPERIMENTS_HOME/software/art/art_illumina must have been built
# 5. $QSIM_EXPERIMENTS_HOME/software/mason/mason must have been built

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
include ../wgsim.mk
include ../art.mk
include ../bowtie2ts.mk

NREADS=4000000

.PHONY: clean
clean:
	rm -f *.fq.gz
	rm -rf *.out

.PHONY: outs
outs: outs_unpaired outs_paired

.PHONY: outs_unpaired
outs_unpaired: \
    r0_bt2s_mason_ill_hg_100.out/DONE \
    r0_bt2s_mason_ill_hg_250.out/DONE \
    r0_bt2s_wgsim_ill_hg_100.out/DONE \
    r0_bt2s_wgsim_ill_hg_250.out/DONE \
    r0_bt2s_art_ill_hg_100.out/DONE \
    r0_bt2s_art_ill_hg_250.out/DONE

.PHONY: outs_paired
outs_paired: \
    r12_bt2s100_mason_ill_hg_100.out/DONE \
    r12_bt2s250_mason_ill_hg_250.out/DONE \
    r12_bt2s100_wgsim_ill_hg_100.out/DONE \
    r12_bt2s250_wgsim_ill_hg_250.out/DONE \
    r12_bt2s100_art_ill_hg_100.out/DONE \
    r12_bt2s250_art_ill_hg_250.out/DONE

.PHONY: reads
reads: \
    r0_mason_ill_hg_100.fq.gz \
    r0_mason_ill_hg_250.fq.gz \
    r0_wgsim_ill_hg_100.fq.gz \
    r0_wgsim_ill_hg_250.fq.gz \
    r0_art_ill_hg_100.fq.gz \
    r0_art_ill_hg_250.fq.gz \
    r1_mason_ill_hg_100.fq.gz \
    r1_mason_ill_hg_250.fq.gz \
    r1_wgsim_ill_hg_100.fq.gz \
    r1_wgsim_ill_hg_250.fq.gz \
    r1_art_ill_hg_100.fq.gz \
    r1_art_ill_hg_250.fq.gz

SUBSAMPLING_ARGS=--subsampling-series 0.015,0.03 --trials=5 --predict-for-training --keep-intermediates --assess-accuracy

#
# Make bowtie2 rule.  Use --sensitive, hg19, and prefix rules with "bt2s"
#
$(eval $(call bt2ts,bt2s,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),--sensitive,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2s100,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 200 -X 400 --sensitive,$(HG_FASTA),$(HG_INDEX),out))
$(eval $(call bt2ts,bt2s250,$(NREADS),$(NREADS),$(SUBSAMPLING_ARGS),-I 500 -X 1000 --sensitive,$(HG_FASTA),$(HG_INDEX),out))

#
# Make the read rules
#
$(eval $(call mason_ill_unp_reads,ill_hg_100,$(HG_FASTA),100,$(NREADS),7724))
$(eval $(call mason_ill_unp_reads,ill_hg_250,$(HG_FASTA),250,$(NREADS),7724))
$(eval $(call wgsim_ill_unp_reads,ill_hg_100,$(HG_FASTA),100,$(NREADS),7724))
$(eval $(call wgsim_ill_unp_reads,ill_hg_250,$(HG_FASTA),250,$(NREADS),7724))

$(eval $(call mason_ill_pair_reads,ill_hg_100,$(HG_FASTA),100,$(NREADS),300,25,7724))
$(eval $(call mason_ill_pair_reads,ill_hg_250,$(HG_FASTA),250,$(NREADS),750,75,7724))
$(eval $(call wgsim_ill_pair_reads,ill_hg_100,$(HG_FASTA),100,$(NREADS),300,25,7724))
$(eval $(call wgsim_ill_pair_reads,ill_hg_250,$(HG_FASTA),250,$(NREADS),750,75,7724))

# For Art, we have to convert desired # reads to a fold coverage.
# Note: for paired-end, we're specifying *fragment* coverage.

# >>> 4000000 * 100 / 3.1e9
# 0.12903225806451613
$(eval $(call art_ill_unp_reads,ill_hg_100,$(HG_FASTA),100,0.129,7724))
# >>> 4000000 * 250 / 3.1e9
# 0.3225806451612903
$(eval $(call art_ill_unp_reads,ill_hg_250,$(HG_FASTA),250,0.323,7724))

# >>> 4000000 * 300 / 3.1e9
# 0.3870967741935484
$(eval $(call art_ill_pair_reads,ill_hg_100,$(HG_FASTA),100,0.387,300,25,7724))
# >>> 4000000 * 750 / 3.1e9
# 0.967741935483871
$(eval $(call art_ill_pair_reads,ill_hg_250,$(HG_FASTA),250,0.968,750,75,7724))
