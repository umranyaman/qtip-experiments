# Must already have a Bowtie 2 index of the FASTA reference somewhere

QSIM=python $(QSIM_HOME)/src/qsim
BT2_QSIM_ARGS=--write-orig-mapq

define bt2ts

r0_$1_%.$8/DONE: r0_%.fq.gz
	mkdir -p $$(basename $$@).temp
	$$(QSIM) --ref $6 \
	       --bt2-exe $$(BOWTIE2) \
	       --index $7 \
	       --sim-unp-min $2 \
	       $$(BT2_TS_ARGS) $$(TS_ARGS) $4 \
	       --output-directory $$(basename $$@) \
	       --temp-directory $$(basename $$@).temp \
	       --U $$< \
	       -- $5 $$(BT2_ARGS) -p 8
	-$$(BOWTIE2) --version > $$(basename $$@)/bt2_version
	-$$(QSIM) --version > $$(basename $$@)/qsim_version
	touch $$@

r12_$1_%.$8/DONE: r1_%.fq.gz
	mkdir -p $$(basename $$@).temp
	$$(QSIM) --ref $6 \
	       --bt2-exe $$(BOWTIE2) \
	       --index $7 \
	       --sim-conc-min $2 --sim-disc-min $3 --sim-bad-end-min $3 \
	       $$(BT2_TS_ARGS) $$(TS_ARGS) $4 \
	       --output-directory $$(basename $$@) \
	       --temp-directory $$(basename $$@).temp \
	       --m1 $$< --m2 $$(<:r1_%=r2_%) \
	       -- $5 $$(BT2_ARGS) -p 8
	-$$(BOWTIE2) --version > $$(basename $$@)/bt2_version
	-$$(QSIM) --version > $$(basename $$@)/qsim_version
	touch $$@

endef
