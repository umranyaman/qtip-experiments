# Must already have a SNAP index of the FASTA reference somewhere

# Arguments:
# 1. Experiment name
# 2. Minimum # unpaired reads to simulate
# 3. Minimum # discordant and bad-end reads to simulate
# 4. Arguments to qtip script
# 5. Arguments for SNAP (in either single or paired mode)
# 6. Arguments for SNAP in single mode
# 7. Arguments for SNAP in paired mode
# 8. Reference to use
# 9. Index to use

QTIP=python $(QTIP_HOME)/src/qtip
SNAP_ARGS+=-=
SNAP_QTIP_ARGS=--write-orig-mapq --write-precise-mapq

define snapts

r0_$1_%.$8/DONE: r0_%.fq.gz
	mkdir -p $$(shell dirname $$(@)).temp
	$$(QTIP) --ref $6 \
	       --snap-exe $$(SNAP) --aligner snap \
	       --index $7.snap \
	       $$(SNAP_QTIP_ARGS) $$(TS_ARGS) $2 \
	       --output-directory $$(shell dirname $$(@)) \
	       --temp-directory $$(shell dirname $$(@)).temp \
	       --U $$< \
	       -- $3 $$(SNAP_ARGS) -t $9 -- $4 -- $5
	$$(SNAP) 2> $$(shell dirname $$(@))/snap_version || true
	$$(QTIP) --version > $$(shell dirname $$(@))/qtip_version || true
	touch $$(@)

r12_$1_%.$8/DONE: r1_%.fq.gz
	mkdir -p $$(shell dirname $$(@)).temp
	$$(QTIP) --ref $6 \
	       --snap-exe $$(SNAP) --aligner snap \
	       --index $7.snap \
	       $$(SNAP_QTIP_ARGS) $$(TS_ARGS) $2 \
	       --output-directory $$(shell dirname $$(@)) \
	       --temp-directory $$(shell dirname $$(@)).temp \
	       --m1 $$< --m2 $$(<:r1_%=r2_%) \
	       -- $3 $$(SNAP_ARGS) -t $9 -- $4 -- $5
	$$(SNAP) 2> $$(shell dirname $$(@))/snap_version || true
	$$(QTIP) --version > $$(shell dirname $$(@))/qtip_version || true
	touch $$(@)

endef
