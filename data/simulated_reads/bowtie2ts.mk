TS=python $(TS_HOME)/bin/ts.py
TS_ARGS=--write-test-data --write-training-data --write-test-distances --write-input-sam --compress-output --verbose --write-timing
BT2_EXTRA_ARGS=--mm

define bt2ts

r0_$1_%.ts: r0_%.fq.gz
	$$(TS) --ref $6 --bt2-exe $$(BOWTIE2) --index $7 --sim-fraction 0.01 --sim-unp-min $2 $$(TS_ARGS) $4 --output-directory $$@ --U $$< -- $5 $$(BT2_ARGS) $$(BT2_EXTRA_ARGS)
	$$(BOWTIE2) --version > $$@/bt2_version
	$$(TS) --version > $$@/ts_version

r12_$1_%.ts: r1_%.fq.gz
	$$(TS) --ref $6 --bt2-exe $$(BOWTIE2) --index $7 --sim-fraction 0.01 --sim-conc-min $2 --sim-disc-min $3 --sim-bad-end-min $3 $$(TS_ARGS) $4 --output-directory $$@ --m1 $$< --m2 $$(<:r1_%=r2_%) -- $5 $$(BT2_ARGS) $$(BT2_EXTRA_ARGS)
	$$(BOWTIE2) --version > $$@/bt2_version
	$$(TS) --version > $$@/ts_version

endef
