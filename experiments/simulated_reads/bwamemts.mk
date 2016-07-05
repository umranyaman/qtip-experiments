# Must already have a BWA index of the FASTA reference somewhere

QSIM=python $(QSIM_HOME)/src/qsim
BWA_TS_ARGS=--write-orig-mapq

define bwamemts

r0_$1_%.$6/DONE: r0_%.fq.gz
	mkdir -p $$(shell dirname $$(@)).temp
	$$(QSIM) --ref $4 \
	       --bwa-exe $$(BWA) --aligner=bwa-mem \
	       --index $5 \
	       $$(BWA_TS_ARGS) $$(TS_ARGS) $2 \
	       --output-directory $$(shell dirname $$(@))\
	       --temp-directory $$(shell dirname $$(@)).temp \
	       --U $$< \
	       -- $3 $$(BWA_ARGS) $$(BWA_EXTRA_ARGS) -t $7
	-$$(BWA) > $$(shell dirname $$(@))/bwa_version 2>&1
	-$$(QSIM) --version > $$(shell dirname $$(@))/qsim_version
	touch $$(@)

r12_$1_%.$6/DONE: r1_%.fq.gz
	mkdir -p $$(shell dirname $$(@)).temp
	$$(QSIM) --ref $4 \
	       --bwa-exe $$(BWA) --aligner=bwa-mem \
	       --index $5 \
	       $$(BWA_TS_ARGS) $$(TS_ARGS) $2 \
	       --output-directory $$(shell dirname $$(@)) \
	       --temp-directory $$(shell dirname $$(@)).temp \
	       --m1 $$< --m2 $$(<:r1_%=r2_%) \
	       -- $3 $$(BWA_ARGS) $$(BWA_EXTRA_ARGS) -t $7
	-$$(BWA) > $$(shell dirname $$(@))/bwa_version 2>&1
	-$$(QSIM) --version > $$(shell dirname $$(@))/qsim_version
	touch $$(@)

endef
