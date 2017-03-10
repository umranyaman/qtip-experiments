#!gmake

# Generate simulated FASTQ files for Illumina-like reads, both
# unpaired and paired-end.  We generate several files, each for a
# different read length from 50 up to 500 nt.  Uses the Mason
# simulator.

# Prerequisites:
# 1. QTIP_EXPERIMENTS_HOME (root of repo) must be set
# 2. $QTIP_EXPERIMENTS_HOME/experiments/refs must have relevant references and indexes
# 3. $QTIP_EXPERIMENTS_HOME/software/mason/mason must have been built

HG_FASTA_NAME=hg38.fa
HG_INDEX_BASE=$(HG_FASTA_NAME)

ifeq (,$(QTIP_EXPERIMENTS_HOME))
$(error "Set QTIP_EXPERIMENTS_HOME, where qtip-experiments working copy is located")
endif
TS_REFS=$(QTIP_EXPERIMENTS_HOME)/experiments/refs
TS_INDEXES=$(QTIP_EXPERIMENTS_HOME)/experiments/refs

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

BOWTIE2_HOME=$(QTIP_EXPERIMENTS_HOME)/software/bowtie2
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
	rm -rf *.out

#.PHONY: outs
#outs: outs_bt2 outs_snap outs_bwamem

#outs_bt2: \
#    r0_bt2s_mason_ill_100_50M.out/DONE \
#    r0_bt2s_mason_ill_250_50M.out/DONE \
#    r12_bt2s100_mason_ill_100_50M.out/DONE \
#    r12_bt2s250_mason_ill_250_50M.out/DONE

#outs_snap: \
#    r0_snap_mason_ill_100_50M.out/DONE \
#    r0_snap_mason_ill_250_50M.out/DONE \
#    r12_snap100_mason_ill_100_50M.out/DONE \
#    r12_snap250_mason_ill_250_50M.out/DONE

#outs_bwamem: \
#    r0_bwamem_mason_ill_100_50M.out/DONE \
#    r0_bwamem_mason_ill_250_50M.out/DONE \
#    r12_bwamem100_mason_ill_100_50M.out/DONE \
#    r12_bwamem250_mason_ill_250_50M.out/DONE

# Problem: mason needs so much memory for paired-end reads, it exhausts RAM on langmead-fs1

.PHONY: reads
reads: \
    r0_mason_ill_100_50M.fq.gz \
    r0_mason_ill_250_50M.fq.gz \
    r1_mason_ill_100_50M.fq.gz \
    r1_mason_ill_250_50M.fq.gz

r1_mason_ill_100_50M.fq.gz: r1_mason_ill_100_1_50M.fq.gz \
                            r1_mason_ill_100_2_50M.fq.gz \
                            r1_mason_ill_100_3_50M.fq.gz \
                            r1_mason_ill_100_4_50M.fq.gz \
                            r1_mason_ill_100_5_50M.fq.gz
	gzip -dc $^ | head -n 200000000 | gzip -c > $@

r1_mason_ill_250_50M.fq.gz: r1_mason_ill_250_1_50M.fq.gz \
                            r1_mason_ill_250_2_50M.fq.gz \
                            r1_mason_ill_250_3_50M.fq.gz \
                            r1_mason_ill_250_4_50M.fq.gz \
                            r1_mason_ill_250_5_50M.fq.gz
	gzip -dc $^ | head -n 200000000 | gzip -c > $@

MK_QTIP_ARGS=--trials=10 --predict-for-training --keep-intermediates --assess-accuracy --skip-rewrite

#$(eval $(call bt2ts,bt2s,1,1,,,$(HG_FASTA),$(HG_INDEX)))
#$(eval $(call bt2ts,bt2s100,1,1,,-I 200 -X 400,$(HG_FASTA),$(HG_INDEX)))
#$(eval $(call bt2ts,bt2s250,1,1,,-I 500 -X 1000,$(HG_FASTA),$(HG_INDEX)))

#$(eval $(call bwamemts,bwamem,1,1,,,$(HG_FASTA),$(HG_INDEX)))
#$(eval $(call bwamemts,bwamem100,1,1,,-I 300$(,)25,$(HG_FASTA),$(HG_INDEX)))
#$(eval $(call bwamemts,bwamem250,1,1,,-I 750$(,)62,$(HG_FASTA),$(HG_INDEX)))

#$(eval $(call snapts,snap,1,1,,,,,$(HG_FASTA),$(HG_INDEX)))
#$(eval $(call snapts,snap100,1,1,,,,-s 100 300,$(HG_FASTA),$(HG_INDEX)))
#$(eval $(call snapts,snap250,1,1,,,,-s 250 750,$(HG_FASTA),$(HG_INDEX)))

# Trim down a bit
$(eval $(call mason_ill_unp_reads,ill_100_50M,$(HG_FASTA),100,50500000,7724))
$(eval $(call mason_ill_unp_reads,ill_250_50M,$(HG_FASTA),250,50500000,7724))

# Cat and trim down a bit
$(eval $(call mason_ill_pair_reads,ill_100_1_50M,$(HG_FASTA),100,10100000,300,100,7724))
$(eval $(call mason_ill_pair_reads,ill_100_2_50M,$(HG_FASTA),100,10100000,300,100,7247))
$(eval $(call mason_ill_pair_reads,ill_100_3_50M,$(HG_FASTA),100,10100000,300,100,2334))
$(eval $(call mason_ill_pair_reads,ill_100_4_50M,$(HG_FASTA),100,10100000,300,100,1923))
$(eval $(call mason_ill_pair_reads,ill_100_5_50M,$(HG_FASTA),100,10100000,300,100,9542))

# Cat and trim down a bit
$(eval $(call mason_ill_pair_reads,ill_250_1_50M,$(HG_FASTA),250,10100000,750,250,7724))
$(eval $(call mason_ill_pair_reads,ill_250_2_50M,$(HG_FASTA),250,10100000,750,250,7247))
$(eval $(call mason_ill_pair_reads,ill_250_3_50M,$(HG_FASTA),250,10100000,750,250,2334))
$(eval $(call mason_ill_pair_reads,ill_250_4_50M,$(HG_FASTA),250,10100000,750,250,1923))
$(eval $(call mason_ill_pair_reads,ill_250_5_50M,$(HG_FASTA),250,10100000,750,250,9542))
