TS=python $(TS_HOME)/bin/ts.py
TS_ARGS=--write-test-data --write-training-data --write-test-distances --write-input-sam --compress-output --verbose

define bwamemts

r0_$1_%.ts: r0_%.fq.gz
	$$(TS) --ref $5 --bwa-exe $$(BWA) --aligner=bwa-mem  --index $6 --num-reads $2 $$(TS_ARGS) $3 --output-directory $$@ --U $$< -- $4 $$(BWA_ARGS) $$(BWA_EXTRA_ARGS)

r12_$1_%.ts: r1_%.fq.gz
	$$(TS) --ref $5 --bwa-exe $$(BWA) --aligner=bwa-mem  --index $6 --num-reads $2 $$(TS_ARGS) $3 --output-directory $$@ --m1 $$< --m2 $$(<:r1_%=r2_%) -- $4 $$(BWA_ARGS) $$(BWA_EXTRA_ARGS) -x $6

endef
