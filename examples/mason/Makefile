#!/bin/sh

TS_PARAMS=--ref $(REF) --fastq --scoring "1,2,6,1,5,3,5,3" --num-reads 10000
GAP=30
BT2_HOME=~/git/bowtie2/bowtie2
BT2_EXE=$(BT2_HOME)/bowtie2
REF=~/bowtie2_indexes/hg19.fa
FA=~/bowtie2_indexes/hg19.fa

all: errs

.PHONY: clean
clean:
	rm -f *.sac *.sat *.sam *.roc *.err *.training

%.sac %.roc: %.sam
	perl ../../bin/correct.pl --gap=$(GAP) -- $< $(<:%.sam=%.roc) $@

%.sat: %.sac
	perl ../../bin/tabulate_sam.pl $< > $@

%.err: %.sat
	Rscript ../../bin/sat2err.R $< | tee $@

define ill_reads

# Tabulated alignments
sats: \
    r0_$1.bt2_vs.sat r0_$1.bt2_s.sat  r0_$1.bt2_f.sat  r0_$1.bt2_vf.sat

# ROC curve
rocs: \
    r0_$1.bt2_vs.roc r0_$1.bt2_s.roc  r0_$1.bt2_f.roc  r0_$1.bt2_vf.roc

# Ranking error
errs: \
    r0_$1.bt2_vs.err r0_$1.bt2_s.err  r0_$1.bt2_f.err  r0_$1.bt2_vf.err

# Alignments at various levels of sensitivity

r0_$1.bt2_vs.sam: r0_$1.fq $$(FA).1.bt2
	python ../../bin/ts.py $$(TS_PARAMS) --save-training r0_$1.bt2_vs.training --U $$< --S $$@ --bt2-exe $$(BT2_EXE) --bt2-args "$$(BT2_ARGS) -x $$(FA) --very-sensitive"

r0_$1.bt2_s.sam:  r0_$1.fq $$(FA).1.bt2
	python ../../bin/ts.py $$(TS_PARAMS) --save-training r0_$1.bt2_s.training  --U $$< --S $$@ --bt2-exe $$(BT2_EXE) --bt2-args "$$(BT2_ARGS) -x $$(FA) --sensitive"

r0_$1.bt2_f.sam:  r0_$1.fq $$(FA).1.bt2
	python ../../bin/ts.py $$(TS_PARAMS) --save-training r0_$1.bt2_f.training  --U $$< --S $$@ --bt2-exe $$(BT2_EXE) --bt2-args "$$(BT2_ARGS) -x $$(FA) --fast"

r0_$1.bt2_vf.sam: r0_$1.fq $$(FA).1.bt2
	python ../../bin/ts.py $$(TS_PARAMS) --save-training r0_$1.bt2_vf.training --U $$< --S $$@ --bt2-exe $$(BT2_EXE) --bt2-args "$$(BT2_ARGS) -x $$(FA) --very-fast"

endef


define ill_reads_ts

# Tabulated alignments
sats: r0_$1.bt2.sat

# ROC curve
rocs: r0_$1.bt2.roc

# Ranking error
errs: r0_$1.bt2.err

r0_$1.bt2.sam: r0_ill_mason_100_50k.fq $$(FA).1.bt2
	python ../../bin/ts.py $4 $$(TS_PARAMS) --save-training r0_$1.bt2.training --U $$< --S $$@ --bt2-exe $$(BT2_EXE) --bt2-args "$$(BT2_ARGS) -x $$(FA)"

endef

#$(eval $(call ill_reads,ill_mason_50_50k,50,50000))
#$(eval $(call ill_reads,ill_mason_75_50k,75,50000))
#$(eval $(call ill_reads,ill_mason_100_50k,100,50000))
#$(eval $(call ill_reads,ill_mason_125_50k,125,50000))
#$(eval $(call ill_reads,ill_mason_150_50k,150,50000))
#$(eval $(call ill_reads,ill_mason_175_50k,175,50000))
#$(eval $(call ill_reads,ill_mason_200_50k,200,50000))
#$(eval $(call ill_reads,ill_mason_225_50k,225,50000))
#$(eval $(call ill_reads,ill_mason_250_50k,250,50000))
#$(eval $(call ill_reads,ill_mason_275_50k,275,50000))
#$(eval $(call ill_reads,ill_mason_300_50k,300,50000))

$(eval $(call ill_reads_ts,ill_mason_075_030_50k,100,50000,--diff-scale=0.75 --num-neighbors=30))
$(eval $(call ill_reads_ts,ill_mason_100_030_50k,100,50000,--diff-scale=1.00 --num-neighbors=30))
$(eval $(call ill_reads_ts,ill_mason_200_030_50k,100,50000,--diff-scale=2.00 --num-neighbors=30))
$(eval $(call ill_reads_ts,ill_mason_500_030_50k,100,50000,--diff-scale=5.00 --num-neighbors=30))

$(eval $(call ill_reads_ts,ill_mason_075_065_50k,100,50000,--diff-scale=0.75 --num-neighbors=65))
$(eval $(call ill_reads_ts,ill_mason_100_065_50k,100,50000,--diff-scale=1.00 --num-neighbors=65))
$(eval $(call ill_reads_ts,ill_mason_200_065_50k,100,50000,--diff-scale=2.00 --num-neighbors=65))
$(eval $(call ill_reads_ts,ill_mason_500_065_50k,100,50000,--diff-scale=5.00 --num-neighbors=65))

$(eval $(call ill_reads_ts,ill_mason_075_100_50k,100,50000,--diff-scale=0.75 --num-neighbors=100))
$(eval $(call ill_reads_ts,ill_mason_100_100_50k,100,50000,--diff-scale=1.00 --num-neighbors=100))
$(eval $(call ill_reads_ts,ill_mason_200_100_50k,100,50000,--diff-scale=2.00 --num-neighbors=100))
$(eval $(call ill_reads_ts,ill_mason_500_100_50k,100,50000,--diff-scale=5.00 --num-neighbors=100))
