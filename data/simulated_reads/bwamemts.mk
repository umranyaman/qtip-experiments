# Must already have a BWA index of the FASTA reference somewhere

TS=python $(TS_HOME)/src/ts.py
TS_ARGS=--compress-output --verbose --write-all --input-reads-simulated

define bwamemts

r0_$1_%.ts: r0_%.fq.gz
	$$(TS) --ref $6 --bwa-exe $$(BWA) --aligner=bwa-mem  --index $7 --sim-fraction 0.01 --sim-unp-min $2 $$(TS_ARGS) $4 --output-directory $$@ --U $$< -- $5 $$(BWA_ARGS) $$(BWA_EXTRA_ARGS)
	-$$(BWA) > $$@/bwa_version 2>&1
	$$(TS) --version > $$@/ts_version

r12_$1_%.ts: r1_%.fq.gz
	$$(TS) --ref $6 --bwa-exe $$(BWA) --aligner=bwa-mem  --index $7 --sim-fraction 0.01 --sim-conc-min $2 --sim-disc-min $3 --sim-bad-end-min $3 $$(TS_ARGS) $4 --output-directory $$@ --m1 $$< --m2 $$(<:r1_%=r2_%) -- $5 $$(BWA_ARGS) $$(BWA_EXTRA_ARGS)
	-$$(BWA) > $$@/bwa_version 2>&1
	$$(TS) --version > $$@/ts_version

endef
