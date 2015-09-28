# Must already have a SNAP index of the FASTA reference somewhere

# Arguments:
# 1. Experiment name
# 2. Minimum # unpaired reads to simulate
# 3. Minimum # discordant and bad-end reads to simulate
# 4. Arguments to ts.py script
# 5. Arguments for SNAP (in either single or paired mode)
# 6. Arguments for SNAP in single mode
# 7. Arguments for SNAP in paired mode
# 8. Reference to use
# 9. Index to use

TS=python $(TS_HOME)/src/ts.py
SNAP_ARGS+=-t 1

define snapts

r0_$1_%.ts: r0_%.fq.gz
	$$(TS) --ref $8 --snap-exe $$(SNAP) --aligner snap --index $9.snap --sim-fraction 0.01 --sim-unp-min $2 $$(SNAP_TS_ARGS) $4 --output-directory $$@ --U $$< -- $5 $$(SNAP_ARGS) -- $6 -- $7
	$$(SNAP) 2> $$@/snap_version
	$$(TS) --version > $$@/ts_version

r12_$1_%.ts: r1_%.fq.gz
	$$(TS) --ref $8 --snap-exe $$(SNAP) --aligner snap --index $9.snap --sim-fraction 0.01 --sim-conc-min $2 --sim-disc-min $3 --sim-bad-end-min $3 $$(SNAP_TS_ARGS) $4 --output-directory $$@ --m1 $$< --m2 $$(<:r1_%=r2_%) -- $5 $$(SNAP_ARGS) -- $6 -- $7
	$$(SNAP) 2> $$@/snap_version
	$$(TS) --version > $$@/ts_version

endef
