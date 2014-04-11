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

.PHONY: all
all: tss

PYTHON=python

BOWTIE2_HOME=$(TS_HOME)/software/bowtie2
BOWTIE2=$(BOWTIE2_HOME)/bowtie2
BT2_ARGS=

$(BOWTIE2):
	$(MAKE) -C $(BOWTIE2_HOME)

include ../../simulated_reads/mason.mk
include ../../simulated_reads/bowtie2ts.mk
include ../../simulated_reads/predict.mk

.PHONY: clean
clean:
	rm -f *.fq.gz
	rm -rf *.ts *.pred

.PHONY: preds
preds: \
    r0_bt2s_chr14_ill_101.pred \
    r0_bt2sl_chr14_ill_101.pred \
    r12_bt2s_chr14_ill_101.pred \
    r12_bt2sl_chr14_ill_101.pred

.PHONY: tss
tss: \
    r0_bt2s_chr14_ill_101.ts \
    r0_bt2sl_chr14_ill_101.ts \
    r12_bt2s_chr14_ill_101.ts \
    r12_bt2sl_chr14_ill_101.ts

.PHONY: reads
reads: \
    r0_chr14_ill_101.fq.gz \
    r1_chr14_ill_101.fq.gz

# Downsample

r0_chr14_ill_101.fq.gz: frag_1.fastq.gz
	$(PYTHON) $(TS_HOME)/bin/fastq_sample.py --in $< --n 1000000 | gzip -c > $@

r1_chr14_ill_101.fq.gz: frag_1.fastq.gz frag_2.fastq.gz
	$(PYTHON) $(TS_HOME)/bin/fastq_sample.py --in $< --n 1000000 | gzip -c > $@
	$(PYTHON) $(TS_HOME)/bin/fastq_sample.py --in $(<:%_1.fastq.gz=%_2.fastq.gz) --n 1000000 | gzip -c > $(@:r1_%=r2_%)

frag_1.fastq.gz:
	wget http://gage.cbcb.umd.edu/data/Hg_chr14/Data.original/frag_1.fastq.gz

frag_2.fastq.gz:
	wget http://gage.cbcb.umd.edu/data/Hg_chr14/Data.original/frag_2.fastq.gz

$(eval $(call bt2ts,bt2s,100000,100000,--correct-chromosomes=chr14,--sensitive,$(HG_FASTA),$(HG_INDEX)))
$(eval $(call bt2ts,bt2sl,100000,100000,--correct-chromosomes=chr14,--sensitive --local,$(HG_FASTA),$(HG_INDEX)))
