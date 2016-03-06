# Must already have a Bowtie 2 index of the FASTA reference somewhere

TS=python $(TS_HOME)/src/qsim
BT2_TS_ARGS=--write-orig-mapq

define bt2ts

r0_$1_%.out: r0_%.fq.gz
	$$(TS) --ref $6 \
	       --bt2-exe $$(BOWTIE2) \
	       --index $7 \
	       --sim-unp-min $2 \
	       $$(BT2_TS_ARGS) $$(TS_ARGS) $4 \
	       --output-directory $$@ \
	       --U $$< \
	       -- $5 $$(BT2_ARGS) -p 8
	-$$(BOWTIE2) --version > $$@/bt2_version
	-$$(TS) --version > $$@/ts_version

r12_$1_%.out: r1_%.fq.gz
	$$(TS) --ref $6 \
	       --bt2-exe $$(BOWTIE2) \
	       --index $7 \
	       --sim-conc-min $2 --sim-disc-min $3 --sim-bad-end-min $3 \
	       $$(BT2_TS_ARGS) $$(TS_ARGS) $4 \
	       --output-directory $$@ \
	       --m1 $$< --m2 $$(<:r1_%=r2_%) \
	       -- $5 $$(BT2_ARGS) -p 8
	-$$(BOWTIE2) --version > $$@/bt2_version
	-$$(TS) --version > $$@/ts_version

endef
