define bt2

r0_$1_%.sam.gz: r0_%.fq.gz
	$$(BOWTIE2) $2 -x $3 -U $$< $$(BT2_ARGS) $$(BT2_EXTRA_ARGS) | gzip -c > $$@

# Look like GNU make doesn't seem to like pattern rules where more than
# one prerequisite uses the pattern.  Otherwise, I'd include r1_%.fq.gz
# as a prerequisite.
r12_$1_%.sam.gz: r1_%.fq.gz
	$$(BOWTIE2) $2 -x $3 -1 $$< -2 $$(<:r1_%=r2_%) $$(BT2_ARGS) $$(BT2_EXTRA_ARGS) | gzip -c > $$@

endef
