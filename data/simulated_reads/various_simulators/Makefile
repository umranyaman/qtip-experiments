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
# 5. $TS_HOME/software/wgsim/wgsim must have been built
# 6. $TS_HOME/software/art/art_illumina must have been built

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
all: outs

BOWTIE2_HOME=$(TS_HOME)/software/bowtie2
BOWTIE2=$(BOWTIE2_HOME)/bowtie2
BT2_ARGS=

$(BOWTIE2):
	$(MAKE) -C $(BOWTIE2_HOME)

include ../mason.mk
include ../wgsim.mk
include ../art.mk
include ../bowtie2ts.mk

.PHONY: clean
clean:
	rm -f *.fq.gz
	rm -rf *.out

.PHONY: outs
outs: outs_unpaired outs_paired

.PHONY: outs_unpaired
outs_unpaired: \
    r0_bt2s_mason_ill_hg_100_1M.out \
    r0_bt2s_mason_ill_hg_250_1M.out \
    r0_bt2s_wgsim_ill_hg_100_1M.out \
    r0_bt2s_wgsim_ill_hg_250_1M.out \
    r0_bt2s_art_ill_hg_100_1M.out \
    r0_bt2s_art_ill_hg_250_1M.out

.PHONY: outs_paired
outs_paired: \
    r12_bt2s100_mason_ill_hg_100_1M.out \
    r12_bt2s250_mason_ill_hg_250_1M.out \
    r12_bt2s100_wgsim_ill_hg_100_1M.out \
    r12_bt2s250_wgsim_ill_hg_250_1M.out \
    r12_bt2s100_art_ill_hg_100_1M.out \
    r12_bt2s250_art_ill_hg_250_1M.out

.PHONY: reads
reads: \
    r0_mason_ill_hg_100_1M.fq.gz \
    r0_mason_ill_hg_250_1M.fq.gz \
    r0_wgsim_ill_hg_100_1M.fq.gz \
    r0_wgsim_ill_hg_250_1M.fq.gz \
    r0_art_ill_hg_100_1M.fq.gz \
    r0_art_ill_hg_250_1M.fq.gz \
    r1_mason_ill_hg_100_1M.fq.gz \
    r1_mason_ill_hg_250_1M.fq.gz \
    r1_wgsim_ill_hg_100_1M.fq.gz \
    r1_wgsim_ill_hg_250_1M.fq.gz \
    r1_art_ill_hg_100_1M.fq.gz \
    r1_art_ill_hg_250_1M.fq.gz

SUBSAMPLING_ARGS=--subsampling-series 0.05,0.1,0.2,0.3,0.4,0.5,1.0 --trials=5 --predict-for-training --keep-intermediates --assess-accuracy

#
# Make bowtie2 rule.  Use --sensitive, hg19, and prefix rules with "bt2s"
#
$(eval $(call bt2ts,bt2s,1000000,1000000,$(SUBSAMPLING_ARGS),--sensitive,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2s100,1000000,1000000,$(SUBSAMPLING_ARGS),-I 200 -X 400 --sensitive,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2s250,1000000,1000000,$(SUBSAMPLING_ARGS),-I 500 -X 1000 --sensitive,$(HG_FASTA),$(HG_INDEX)))

#
# Make the read rules
#
$(eval $(call mason_ill_unp_reads,ill_hg_100_1M,$(HG_FASTA),100,1000000,7724))
$(eval $(call mason_ill_unp_reads,ill_hg_250_1M,$(HG_FASTA),250,1000000,7724))
$(eval $(call wgsim_ill_unp_reads,ill_hg_100_1M,$(HG_FASTA),100,1000000,7724))
$(eval $(call wgsim_ill_unp_reads,ill_hg_250_1M,$(HG_FASTA),250,1000000,7724))

$(eval $(call mason_ill_pair_reads,ill_hg_100_1M,$(HG_FASTA),100,1000000,300,25,7724))
$(eval $(call mason_ill_pair_reads,ill_hg_250_1M,$(HG_FASTA),250,1000000,750,75,7724))
$(eval $(call wgsim_ill_pair_reads,ill_hg_100_1M,$(HG_FASTA),100,1000000,300,25,7724))
$(eval $(call wgsim_ill_pair_reads,ill_hg_250_1M,$(HG_FASTA),250,1000000,750,75,7724))

# For Art, we have to convert desired # reads to a fold coverage.
# Note: for paired-end, we're specifying *fragment* coverage.

# >>> 1000000 * 100 / 3.1e9
# 0.032258064516129032
$(eval $(call art_ill_unp_reads,ill_hg_100_1M,$(HG_FASTA),100,0.033,7724))
# >>> 1000000 * 250 / 3.1e9
# 0.08064516129032258
$(eval $(call art_ill_unp_reads,ill_hg_250_1M,$(HG_FASTA),250,0.081,7724))

# >>> 1000000 * 300 / 3.1e9
# 0.0967741935483871
$(eval $(call art_ill_pair_reads,ill_hg_100_1M,$(HG_FASTA),100,0.1,300,25,7724))
# >>> 1000000 * 750 / 3.1e9
# 0.24193548387096774
$(eval $(call art_ill_pair_reads,ill_hg_250_1M,$(HG_FASTA),250,0.25,750,75,7724))
