# Must already have a BWA index of the FASTA reference somewhere

TS=python $(TS_HOME)/src/qsim
BWA_TS_ARGS=--write-orig-mapq

define bwamemts

r0_$1_%.out: r0_%.fq.gz
	$$(TS) --ref $6 \
	       --bwa-exe $$(BWA) --aligner=bwa-mem \
	       --index $7 \
	       --sim-unp-min $2 \
	       $$(BWA_TS_ARGS) $$(TS_ARGS) $4 \
	       --output-directory $$@ \
	       --U $$< \
	       -- $5 $$(BWA_ARGS) $$(BWA_EXTRA_ARGS) -t 8
	-$$(BWA) > $$@/bwa_version 2>&1
	-$$(TS) --version > $$@/ts_version

r12_$1_%.out: r1_%.fq.gz
	$$(TS) --ref $6 \
	       --bwa-exe $$(BWA) --aligner=bwa-mem \
	       --index $7 \
	       --sim-conc-min $2 --sim-disc-min $3 --sim-bad-end-min $3 \
	       $$(BWA_TS_ARGS) $$(TS_ARGS) $4 \
	       --output-directory $$@ \
	       --m1 $$< --m2 $$(<:r1_%=r2_%) \
	       -- $5 $$(BWA_ARGS) $$(BWA_EXTRA_ARGS) -t 8
	-$$(BWA) > $$@/bwa_version 2>&1
	-$$(TS) --version > $$@/ts_version

endef
