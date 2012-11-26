#!gmake

##
# Makefile
# 
# Starting from scratch, obtains wgsim, source for Bowtie 2 and BWA, builds all
# three, downloads reference, indexes it with both tools' indexers, simulates
# reads with wgsim, aligns them, and generates ROCs and R-friendly fixed-witth
# SAM tables.
#

SVN=$(shell which svn)
MASON=seqan-trunk/build/Release/core/apps/mason/mason
TOOLBOX=toolbox

#
# Aligners
#

# Bowtie 2
BT2_VER_SM=2.0.0-beta7
BT2_VER=bowtie2-$(BT2_VER_SM)
BT2_AR=$(BT2_VER)-source.zip
BT2_URL=http://sourceforge.net/projects/bowtie-bio/files/bowtie2/$(BT2_VER_SM)/$(BT2_AR)/download
BT2_ARGS=-t --mapq-extra -I 100 -X 550 --read-times

# Bowtie
BT_VER_SM=0.12.8
BT_VER=bowtie-$(BT_VER_SM)
BT_AR=$(BT_VER)-src.zip
BT_URL=http://sourceforge.net/projects/bowtie-bio/files/bowtie/$(BT_VER_SM)/$(BT_AR)/download
BT_ARGS=-t -I 100 -X 550 --sam

# BWA
BWA_VER=bwa-0.5.9
BWA_AR=$(BWA_VER).tar.bz2
BWA_URL=http://sourceforge.net/projects/bio-bwa/files/$(BWA_AR)/download

# SNAP
SNAP_VER=

#
# Other tools
#

# # CPUs available to subtasks, like the simulation scripts
NCPUS=1

# Reference genome FASTA file
FA=hg19.fa

# RepeatMasker file
REP=hg19.fa.out.gz
REP_URL=http://www.repeatmasker.org/genomes/hg19/RepeatMasker-rm330-db20120124

# Repeat annotation script
REP_SCR=$(TOOLBOX)/python/bin/rep.py

%.sat %.roc: %.sam $(REP)
	perl correct.pl $< $(<:%.sam=%.roc) .tmp.$<
	perl tabulate_sam.pl .tmp.$< | python $(REP_SCR) --rep=$(REP) --sat-in=/dev/stdin --sat-out=$@
	rm -f .tmp.$<

.PHONY: reads

define ill_reads

ill_reads: r0_$1.fq r1_$1.fq

reads: r0_$1.fq r1_$1.fq

alignments: \
    r0_$1.bt_v.sat r0_$1.bt_default.sat r0_$1.bt_l28n2e100.sat r0_$1.bt_l28n2e150.sat r0_$1.bt_l28n2e200.sat r0_$1.bt_l28n2e250.sat \
    r0_$1.bt2_vs.sat r0_$1.bt2_s.sat r0_$1.bt2_f.sat r0_$1.bt2_vf.sat \
    r0_$1.bwa.sat \
    r12_$1.bt_default.sat r12_$1.bt_l28n2e250.sat \
    r12_$1.bt2_vs.sat r12_$1.bt2_s.sat r12_$1.bt2_f.sat r12_$1.bt2_vf.sat \
    r12_$1.bwa.sat

rocs: \
    r0_$1.bt_v.roc r0_$1.bt_default.roc r0_$1.bt_l28n2e100.roc r0_$1.bt_l28n2e150.roc r0_$1.bt_l28n2e200.roc r0_$1.bt_l28n2e250.roc \
    r0_$1.bt2_vs.roc r0_$1.bt2_s.roc  r0_$1.bt2_f.roc  r0_$1.bt2_vf.roc \
    r0_$1.bwa.roc \
    r12_$1.bt_default.roc r12_$1.bt_l28n2e250.roc \
    r12_$1.bt2_vs.roc r12_$1.bt2_s.roc r12_$1.bt2_f.roc r12_$1.bt2_vf.roc \
    r12_$1.bwa.roc

r0_$1.bt_v.sam: r0_$1.fq $$(BT_VER)/bowtie $$(FA).1.ebwt
	$$(BT_VER)/bowtie -v 2 -M 1 --best $$(BT_ARGS) $$(FA) $$< > $$@

r0_$1.bt_l28n2e100.sam: r0_$1.fq $$(BT_VER)/bowtie $$(FA).1.ebwt
	$$(BT_VER)/bowtie -l 28 -n 2 -e 100 -M 1 --best $$(BT_ARGS) $$(FA) $$< > $$@

r0_$1.bt_default.sam: r0_$1.fq $$(BT_VER)/bowtie $$(FA).1.ebwt
	$$(BT_VER)/bowtie -M 1 --best $$(BT_ARGS) $$(FA) $$< > $$@

r0_$1.bt_l28n2e150.sam: r0_$1.fq $$(BT_VER)/bowtie $$(FA).1.ebwt
	$$(BT_VER)/bowtie -l 28 -n 2 -e 150 -M 1 --best $$(BT_ARGS) $$(FA) $$< > $$@

r0_$1.bt_l28n2e200.sam: r0_$1.fq $$(BT_VER)/bowtie $$(FA).1.ebwt
	$$(BT_VER)/bowtie -l 28 -n 2 -e 200 -M 1 --best $$(BT_ARGS) $$(FA) $$< > $$@

r0_$1.bt_l28n2e250.sam: r0_$1.fq $$(BT_VER)/bowtie $$(FA).1.ebwt
	$$(BT_VER)/bowtie -l 28 -n 2 -e 250 -M 1 --best $$(BT_ARGS) $$(FA) $$< > $$@

r0_$1.bt2_vs.sam: r0_$1.fq $$(BT2_VER)/bowtie2-align $$(FA).1.bt2
	$$(BT2_VER)/bowtie2 $$(BT2_ARGS) --very-sensitive -x $$(FA) -U $$< --met 1 --met-file $$(@:%.sam=%.met) > $$@

r0_$1.bt2_s.sam: r0_$1.fq $$(BT2_VER)/bowtie2-align $$(FA).1.bt2
	$$(BT2_VER)/bowtie2 $$(BT2_ARGS) --sensitive -x $$(FA) -U $$< --met 1 --met-file $$(@:%.sam=%.met) > $$@

r0_$1.bt2_f.sam: r0_$1.fq $$(BT2_VER)/bowtie2-align $$(FA).1.bt2
	$$(BT2_VER)/bowtie2 $$(BT2_ARGS) --fast -x $$(FA) -U $$< --met 1 --met-file $$(@:%.sam=%.met) > $$@

r0_$1.bt2_vf.sam: r0_$1.fq $$(BT2_VER)/bowtie2-align $$(FA).1.bt2
	$$(BT2_VER)/bowtie2 $$(BT2_ARGS) --very-fast -x $$(FA) -U $$< --met 1 --met-file $$(@:%.sam=%.met) > $$@

r0_$1.bwa.sam: r0_$1.fq $$(BWA_VER)/bwa $$(FA).amb
	$$(BWA_VER)/bwa aln $$(FA) $$< > .$$<.sai
	$$(BWA_VER)/bwa samse $$(FA) .$$<.sai $$< > $$@

r0_$1.bwasw.sam: r0_$1.fq $$(BWA_VER)/bwa $$(FA).amb
	$$(BWA_VER)/bwa bwasw $$(FA) $$< > $$@

r0_$1.fq: $$(MASON) mason.pl $$(FA)
	perl mason.pl --mason-exe=$$(MASON) --num-cpus=$(NCPUS) --num-reads=$3 --out-fq-1=$$@                            --out-sam=$$(@:%.fq=%.sam) --out-tsv=$$(@:%.fq=%.tsv) -- $$(FA) -- illumina -hn 2 -i -s 322 -sq -n $2

r12_$1.bt_default.sam: r1_$1.fq $$(BT_VER)/bowtie $$(FA).1.ebwt
	$$(BT_VER)/bowtie $$(BT_ARGS) -M 1 --best $$(FA) -1 $$< -2 $$(<:r1_%=r2_%) > $$@

r12_$1.bt_l28n2e250.sam: r1_$1.fq $$(BT_VER)/bowtie $$(FA).1.ebwt
	$$(BT_VER)/bowtie $$(BT_ARGS) -l 28 -n 2 -e 250 -M 1 --best $$(FA) -1 $$< -2 $$(<:r1_%=r2_%) > $$@

r12_$1.bt2_vs.sam: r1_$1.fq $$(BT2_VER)/bowtie2-align $$(FA).1.bt2
	$$(BT2_VER)/bowtie2 $$(BT2_ARGS) --very-sensitive -x $$(FA) -1 $$< -2 $$(<:r1_%=r2_%) --met 1 --met-file $$(@:%.sam=%.met) > $$@

r12_$1.bt2_s.sam: r1_$1.fq $$(BT2_VER)/bowtie2-align $$(FA).1.bt2
	$$(BT2_VER)/bowtie2 $$(BT2_ARGS) --sensitive -x $$(FA) -1 $$< -2 $$(<:r1_%=r2_%) --met 1 --met-file $$(@:%.sam=%.met) > $$@

r12_$1.bt2_f.sam: r1_$1.fq $$(BT2_VER)/bowtie2-align $$(FA).1.bt2
	$$(BT2_VER)/bowtie2 $$(BT2_ARGS) --fast -x $$(FA) -1 $$< -2 $$(<:r1_%=r2_%) --met 1 --met-file $$(@:%.sam=%.met) > $$@

r12_$1.bt2_vf.sam: r1_$1.fq $$(BT2_VER)/bowtie2-align $$(FA).1.bt2
	$$(BT2_VER)/bowtie2 $$(BT2_ARGS) --very-fast -x $$(FA) -1 $$< -2 $$(<:r1_%=r2_%) --met 1 --met-file $$(@:%.sam=%.met) > $$@

r12_$1.bwa.sam: r1_$1.fq $$(BWA_VER)/bwa $$(FA).amb
	$$(BWA_VER)/bwa aln $$(FA) $$< > .$$<.sai
	$$(BWA_VER)/bwa aln $$(FA) $$(<:r1_%=r2_%) > .$$(<:r1_%=r2_%).sai
	$$(BWA_VER)/bwa sampe $$(FA) .$$<.sai .$$(<:r1_%=r2_%).sai $$< $$(<:r1_%=r2_%) > $$@

r1_$1.fq: $(MASON) mason.pl $(FA)
	perl mason.pl --mason-exe=$$(MASON) --num-cpus=$(NCPUS) --num-reads=$3 --out-fq-1=$$@ --out-fq-2=$$(@:r1_%=r2_%) --out-sam=$$(@:%.fq=%.sam) --out-tsv=$$(@:%.fq=%.tsv) -- $$(FA) -- illumina -hn 2 -i -s 535 -sq -mp -rn 2 -ll 375 -le 100 -n $2

endef

define fff_reads

fff_reads: r0_$1.fq

reads: r0_$1.fq

alignments: r0_$1.bt2l_vs.sat r0_$1.bt2l_s.sat r0_$1.bt2l_s2.sat r0_$1.bt2l_s3.sat r0_$1.bt2l_f.sat r0_$1.bt2l_vf.sat \
            r0_$1.bwa.sat

rocs: r0_$1.bt2l_vs.roc r0_$1.bt2l_s.roc r0_$1.bt2l_s2.roc r0_$1.bt2l_s3.roc r0_$1.bt2l_f.roc r0_$1.bt2l_vf.roc \
      r0_$1.bwa.roc

r0_$1.bt2l_vs.sam: r0_$1.fq $$(BT2_VER)/bowtie2-align $$(FA).1.bt2
	$$(BT2_VER)/bowtie2 $$(BT2_ARGS) --very-sensitive --local --bwa-sw-like -x $$(FA) -U $$< --met 1 --met-file $$(@:%.sam=%.met) > $$@

r0_$1.bt2l_s.sam: r0_$1.fq $$(BT2_VER)/bowtie2-align $$(FA).1.bt2
	$$(BT2_VER)/bowtie2 $$(BT2_ARGS) --sensitive --local --bwa-sw-like -x $$(FA) -U $$< --met 1 --met-file $$(@:%.sam=%.met) > $$@

r0_$1.bt2l_f.sam: r0_$1.fq $$(BT2_VER)/bowtie2-align $$(FA).1.bt2
	$$(BT2_VER)/bowtie2 $$(BT2_ARGS) --fast --local --bwa-sw-like -x $$(FA) -U $$< --met 1 --met-file $$(@:%.sam=%.met) > $$@

r0_$1.bt2l_vf.sam: r0_$1.fq $$(BT2_VER)/bowtie2-align $$(FA).1.bt2
	$$(BT2_VER)/bowtie2 $$(BT2_ARGS) --very-fast --local --bwa-sw-like -x $$(FA) -U $$< --met 1 --met-file $$(@:%.sam=%.met) > $$@

r0_$1.bt2l_s2.sam: r0_$1.fq $$(BT2_VER)/bowtie2-align $$(FA).1.bt2
	$$(BT2_VER)/bowtie2 $$(BT2_ARGS) --local --score-min=G,20,8 -x $$(FA) -U $$< --met 1 --met-file $$(@:%.sam=%.met) > $$@

r0_$1.bt2l_s3.sam: r0_$1.fq $$(BT2_VER)/bowtie2-align $$(FA).1.bt2
	$$(BT2_VER)/bowtie2 $$(BT2_ARGS) --local -x $$(FA) -U $$< --met 1 --met-file $$(@:%.sam=%.met) > $$@

r0_$1.bwa.sam: r0_$1.fq $$(BWA_VER)/bwa $$(FA).amb
	$$(BWA_VER)/bwa aln $$(FA) $$< > .$$<.sai
	$$(BWA_VER)/bwa samse $$(FA) .$$<.sai $$< > $$@

r0_$1.bwasw.sam: r0_$1.fq $$(BWA_VER)/bwa $$(FA).amb
	$$(BWA_VER)/bwa bwasw $$(FA) $$< > $$@

r0_$1.fq: $$(MASON) mason.pl $$(FA)
	perl mason.pl --mason-exe=$$(MASON) --num-cpus=$(NCPUS) --num-reads=$3 --out-fq-1=$$@                            --out-sam=$$(@:%.fq=%.sam) --out-tsv=$$(@:%.fq=%.tsv) -- $$(FA) -- 454 -hn 2 -i -s 322 -sq -nm $2 $4

endef

#
# Mason takes a surprisingly long time to run.  If you want to generate all the
# collections of 100K reads below in less than a day, you'll probably have to
# do the runs for the longer read lengths on many cores, using the mason.pl
# script directly.
#

$(eval $(call ill_reads,ill_75_100k,75,100000))
$(eval $(call ill_reads,ill_100_100k,100,100000))
$(eval $(call ill_reads,ill_125_100k,125,100000))
$(eval $(call ill_reads,ill_150_100k,150,100000))
$(eval $(call ill_reads,ill_175_100k,175,100000))
$(eval $(call fff_reads,fff_250_100k,250,100000))
$(eval $(call fff_reads,fff_400_100k,400,100000))
$(eval $(call fff_reads,fffn_250_100k,250,100000,-k 0.3 -bm 0.4 -bs 0.2))
$(eval $(call fff_reads,fffn_400_100k,400,100000,-k 0.3 -bm 0.4 -bs 0.2))
#$(eval $(call fff_reads,fffn_5000_50,5000,50,-k 0.3 -bm 0.4 -bs 0.2))
#$(eval $(call fff_reads,fffn_10000_50,10000,50,-k 0.3 -bm 0.4 -bs 0.2))
#$(eval $(call fff_reads,fffn_50000_20,50000,20,-k 0.3 -bm 0.4 -bs 0.2))

#
# Mason simulator - assumues cmake & gcc tools are already installed
#

$(MASON): seqan-trunk/build
	cd seqan-trunk/build && mkdir -p Release && cd Release && cmake -DCMAKE_BUILD_TYPE=Release ../.. && $(MAKE) mason

seqan-trunk/build: $(SVN)
	$(SVN) co http://svn.mi.fu-berlin.de/seqan/trunk/seqan seqan-trunk
	mkdir -p $@

#
# Indexes
#

.PHONY: indexes
indexes: $(FA).amb $(FA).1.bt2 $(FA).1.ebwt

$(FA).1.bt2: $(BT2_VER)/bowtie2-build $(FA)
	$< $(FA) $(FA)

$(FA).1.ebwt: $(BT_VER)/bowtie-build $(FA)
	$< $(FA) $(FA)

$(FA).amb: $(BWA_VER)/bwa $(FA)
	$< index -a bwtsw $(FA)

$(FA): get_hg19.sh
	sh get_hg19.sh

$(REP):
	wget $(REP_URL)/$(REP)

#
# Tools
#

.PHONY: tools
tools: $(BT2_VER)/bowtie2-align $(BT2_VER)/bowtie2-build \
       $(BT_VER)/bowtie $(BT_VER)/bowtie-build \
       $(BWA_VER)/bwa $(REP_SCR)

$(BT_VER)/bowtie: $(BT_VER)/Makefile
	$(MAKE) -C $(BT_VER) bowtie

$(BT_VER)/bowtie-build: $(BT_VER)/Makefile
	$(MAKE) -C $(BT_VER) bowtie-build

$(BT_VER)/Makefile:
	wget $(BT_URL)
	unzip -u $(BT_AR)

$(BT2_VER)/bowtie2-align: $(BT2_VER)/Makefile
	$(MAKE) -C $(BT2_VER) bowtie2-align

$(BT2_VER)/bowtie2-build: $(BT2_VER)/Makefile
	$(MAKE) -C $(BT2_VER) bowtie2-build

$(BT2_VER)/Makefile:
	wget $(BT2_URL)
	unzip -u $(BT2_AR)

$(BWA_VER)/bwa: $(BWA_VER)
	make -C $<

$(BWA_VER):
	wget $(BWA_URL)
	tar xvfj $(BWA_AR)

$(REP_SCR):
	sh co_toolbox.sh
