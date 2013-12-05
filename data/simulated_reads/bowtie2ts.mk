TS=pypy $(TS_HOME)/bin/ts.py
TS_ARGS=--write-test-data --write-training-data --compress-output

define bt2ts

r0_$1_%.ts: r0_%.fq.gz
	$$(TS) --ref $4 --bt2-exe $$(BOWTIE2) --num-reads $2 $$(TS_ARGS) --output-directory $$@ --U $$< -- $3 $$(BT2_ARGS) $$(BT2_EXTRA_ARGS) -x $5

r12_$1_%.ts: r1_%.fq.gz
	$$(TS) --ref $4 --bt2-exe $$(BOWTIE2) --num-reads $2 $$(TS_ARGS) --output-directory $$@ --m1 $$< --m2 $$(<:r1_%=r2_%) -- $3 $$(BT2_ARGS) $$(BT2_EXTRA_ARGS) -x $5

endef
