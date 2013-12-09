TS=pypy $(TS_HOME)/bin/ts.py
TS_ARGS=--write-test-data --write-training-data --write-test-distances --write-input-sam --compress-output --verbose

define bt2ts

r0_$1_%.ts: r0_%.fq.gz
	$$(TS) --ref $5 --bt2-exe $$(BOWTIE2) --num-reads $2 $$(TS_ARGS) $3 --output-directory $$@ --U $$< -- $4 $$(BT2_ARGS) $$(BT2_EXTRA_ARGS) -x $6

r12_$1_%.ts: r1_%.fq.gz
	$$(TS) --ref $5 --bt2-exe $$(BOWTIE2) --num-reads $2 $$(TS_ARGS) $3 --output-directory $$@ --m1 $$< --m2 $$(<:r1_%=r2_%) -- $4 $$(BT2_ARGS) $$(BT2_EXTRA_ARGS) -x $6

endef
