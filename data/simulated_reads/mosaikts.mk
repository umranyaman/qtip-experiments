# Must already have a "built" version of the FASTA reference somewhere

TS=python -W error $(TS_HOME)/bin/ts.py
TS_ARGS=--compress-output --verbose --write-all

r0_%.fq.mkb: r0_%.fq.gz
	$(MOSAIK_BUILD) -q $< -out $@ -st illumina

r12_%.fq.mkb: r1_%.fq.gz r2_%.fq.gz
	$(MOSAIK_BUILD) -q $< -q2 $(<:r1_%=r2_%) -out $@ -st illumina

define mosaikts

r0_$1_%.ts: r0_%.fq.mkb
	$$(TS) --ref $6 --mosaik-align-exe $$(MOSAIK_ALIGN) --mosaik-build-exe $$(MOSAIK_BUILD) --aligner=mosaik \
           --index $7.mosaik --sim-fraction 0.01 --sim-unp-min $2 $$(TS_ARGS) $4 --output-directory $$@ --U $$< \
           -- $5 $$(MOSAIK_ALIGN_ARGS) \
                 -annpe $(MOSAIK_HOME)/MOSAIK/networkFile/2.1.26.pe.100.0065.ann \
                 -annse $(MOSAIK_HOME)/MOSAIK/networkFile/2.1.26.se.100.005.ann
	-$$(MOSAIK_ALIGN) > $$@/mosaik_align_version 2>&1
	-$$(MOSAIK_BUILD) > $$@/mosaik_build_version 2>&1
	$$(TS) --version > $$@/ts_version

r12_$1_%.ts: r12_%.fq.mkb
	$$(TS) --ref $6 --bwa-exe $$(BWA) --aligner=bwa-mem  --index $7.mosaik --sim-fraction 0.01 \
	       --sim-conc-min $2 --sim-disc-min $3 --sim-bad-end-min $3 $$(TS_ARGS) $4 \
	       --output-directory $$@ --m1 $$< --m2 $$(<:r1_%=r2_%) \
           -- $5 $$(MOSAIK_ALIGN_ARGS) \
                 -annpe $(MOSAIK_HOME)/MOSAIK/networkFile/2.1.26.pe.100.0065.ann \
                 -annse $(MOSAIK_HOME)/MOSAIK/networkFile/2.1.26.se.100.005.ann
	-$$(MOSAIK_ALIGN) > $$@/mosaik_align_version 2>&1
	-$$(MOSAIK_BUILD) > $$@/mosaik_build_version 2>&1
	$$(TS) --version > $$@/ts_version

endef
