# Must already have a BWA index of the FASTA reference somewhere

QSIM=python $(QSIM_HOME)/src/qsim
BWA_TS_ARGS=--write-orig-mapq

define bwamemts

r0_$1_%.$8/DONE: r0_%.fq.gz
	mkdir -p $$(shell dirname $$(@)).temp
	$$(QSIM) --ref $6 \
	       --bwa-exe $$(BWA) --aligner=bwa-mem \
	       --index $7 \
	       --sim-unp-min $2 \
	       $$(BWA_TS_ARGS) $$(TS_ARGS) $4 \
	       --output-directory $$(shell dirname $$(@))\
	       --temp-directory $$(shell dirname $$(@)).temp \
	       --U $$< \
	       -- $5 $$(BWA_ARGS) $$(BWA_EXTRA_ARGS) -t $9
	-$$(BWA) > $$(shell dirname $$(@))/bwa_version 2>&1
	-$$(QSIM) --version > $$(shell dirname $$(@))/qsim_version
	touch $$(@)

r12_$1_%.$8/DONE: r1_%.fq.gz
	mkdir -p $$(shell dirname $$(@)).temp
	$$(QSIM) --ref $6 \
	       --bwa-exe $$(BWA) --aligner=bwa-mem \
	       --index $7 \
	       --sim-conc-min $2 --sim-disc-min $3 --sim-bad-end-min $3 \
	       $$(BWA_TS_ARGS) $$(TS_ARGS) $4 \
	       --output-directory $$(shell dirname $$(@)) \
	       --temp-directory $$(shell dirname $$(@)).temp \
	       --m1 $$< --m2 $$(<:r1_%=r2_%) \
	       -- $5 $$(BWA_ARGS) $$(BWA_EXTRA_ARGS) -t $9
	-$$(BWA) > $$(shell dirname $$(@))/bwa_version 2>&1
	-$$(QSIM) --version > $$(shell dirname $$(@))/qsim_version
	touch $$(@)

endef
