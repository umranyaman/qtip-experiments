# Must already have a SNAP index of the FASTA reference somewhere

# Arguments:
# 1. Experiment name
# 2. Minimum # unpaired reads to simulate
# 3. Minimum # discordant and bad-end reads to simulate
# 4. Arguments to qsim script
# 5. Arguments for SNAP (in either single or paired mode)
# 6. Arguments for SNAP in single mode
# 7. Arguments for SNAP in paired mode
# 8. Reference to use
# 9. Index to use

QSIM=python $(QSIM_HOME)/src/qsim
SNAP_ARGS+=-=
SNAP_TS_ARGS=--write-orig-mapq

define snapts

r0_$1_%.$(10)/DONE: r0_%.fq.gz
	mkdir -p $$(shell dirname $$(@)).temp
	$$(QSIM) --ref $8 \
	       --snap-exe $$(SNAP) --aligner snap \
	       --index $9.snap \
	       --sim-unp-min $2 \
	       $$(SNAP_TS_ARGS) $$(TS_ARGS) $4 \
	       --output-directory $$(shell dirname $$(@)) \
	       --temp-directory $$(shell dirname $$(@)).temp \
	       --U $$< \
	       -- $5 $$(SNAP_ARGS) -t $(11) -- $6 -- $7
	-$$(SNAP) 2> $$(shell dirname $$(@))/snap_version
	-$$(QSIM) --version > $$(shell dirname $$(@))/qsim_version
	touch $$(@)

r12_$1_%.$(10)/DONE: r1_%.fq.gz
	mkdir -p $$(shell dirname $$(@)).temp
	$$(QSIM) --ref $8 \
	       --snap-exe $$(SNAP) --aligner snap \
	       --index $9.snap \
	       --sim-conc-min $2 --sim-disc-min $3 --sim-bad-end-min $3 \
	       $$(SNAP_TS_ARGS) $$(TS_ARGS) $4 \
	       --output-directory $$(shell dirname $$(@)) \
	       --temp-directory $$(shell dirname $$(@)).temp \
	       --m1 $$< --m2 $$(<:r1_%=r2_%) \
	       -- $5 $$(SNAP_ARGS) -t $(11) -- $6 -- $7
	-$$(SNAP) 2> $$(shell dirname $$(@))/snap_version
	-$$(QSIM) --version > $$(shell dirname $$(@))/qsim_version
	touch $$(@)

endef
