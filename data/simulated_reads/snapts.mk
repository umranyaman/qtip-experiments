# Must already have a SNAP index of the FASTA reference somewhere

TS=python $(TS_HOME)/bin/ts.py
TS_ARGS=--compress-output --verbose --write-all
SNAP_EXTRA_ARGS=

define snapts

r0_$1_%.ts: r0_%.fq.gz
	$$(TS) --ref $6 --snap-exe $$(SNAP) --aligner snap --index $7.snap --sim-fraction 0.01 --sim-unp-min $2 $$(TS_ARGS) $4 --output-directory $$@ --U $$< -- $5 $$(SNAP_ARGS) $$(SNAP_EXTRA_ARGS)
	$$(SNAP) 2> $$@/snap_version
	$$(TS) --version > $$@/ts_version

r12_$1_%.ts: r1_%.fq.gz
	$$(TS) --ref $6 --snap-exe $$(SNAP) --aligner snap --index $7.snap --sim-fraction 0.01 --sim-conc-min $2 --sim-disc-min $3 --sim-bad-end-min $3 $$(TS_ARGS) $4 --output-directory $$@ --m1 $$< --m2 $$(<:r1_%=r2_%) -- $5 $$(SNAP_ARGS) $$(SNAP_EXTRA_ARGS)
	$$(SNAP) 2> $$@/snap_version
	$$(TS) --version > $$@/ts_version

endef
