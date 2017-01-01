# Must already have a Bowtie 2 index of the FASTA reference somewhere

QTIP=python $(QTIP_HOME)/src/qtip
BT2_QTIP_ARGS=--write-orig-mapq --write-precise-mapq

define bt2ts

r0_$1_%.$6/DONE: r0_%.fq.gz
	mkdir -p $$(shell dirname $$(@)).temp
	$$(QTIP) --ref $4 \
	       --bt2-exe $$(BOWTIE2) \
	       --index $5 \
	       $$(BT2_QTIP_ARGS) $$(TS_ARGS) $2 \
	       --output-directory $$(shell dirname $$(@)) \
	       --temp-directory $$(shell dirname $$(@)).temp \
	       --U $$< \
	       -- $3 $$(BT2_ARGS) -p $7
	$$(BOWTIE2) --version > $$(shell dirname $$(@))/bt2_version || true
	$$(QTIP) --version > $$(shell dirname $$(@))/qtip_version || true
	touch $$(@)

r12_$1_%.$6/DONE: r1_%.fq.gz
	mkdir -p $$(shell dirname $$(@)).temp
	$$(QTIP) --ref $4 \
	       --bt2-exe $$(BOWTIE2) \
	       --index $5 \
	       $$(BT2_QTIP_ARGS) $$(TS_ARGS) $2 \
	       --output-directory $$(shell dirname $$(@)) \
	       --temp-directory $$(shell dirname $$(@)).temp \
	       --m1 $$< --m2 $$(<:r1_%=r2_%) \
	       -- $3 $$(BT2_ARGS) -p $7
	$$(BOWTIE2) --version > $$(shell dirname $$(@))/bt2_version || true
	$$(QTIP) --version > $$(shell dirname $$(@))/qtip_version || true
	touch $$(@)

endef
