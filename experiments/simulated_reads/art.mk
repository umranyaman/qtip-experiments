# art_illumina [options] -i <DNA_reference_file> -l <read_length> -f <fold_coverage> -o <outFile_prefix>
#
# Relevant options:
#  -i   --in       the filename of input DNA/RNA reference
#  -mp  --matepair indicate a mate-pair read simulation
#  -p   --paired   indicate a paired-end read simulation
#       NOTE: art will automatically switch to a mate-pair read simulation if the given mean fragment size >= 2000
#  -l   --len      the length of reads to be simulated
#  -f   --fcov     the fold of read coverage to be simulated
#  -m   --mflen    the mean size of DNA fragments for paired-end simulations
#  -s   --sdev     the standard deviation of DNA fragment size for paired-end simulations.
#  -na  --noALN    do not output ALN alignment file
#  -rs  --rndSend  the seed for random number generator (default: system time in second)

# Use art to generate unpaired Illumina-like reads.  Macro takes
# parameters: (1) batch name, (2) reference genome, (3) read length in
# batch, (4) target # reads, (5) fold coverage, (6) pseudo-random seed.
define art_ill_unp_reads

# Generate unpaired reads
r0_art_$1.fq.gz: $$(FA) $$(QTIP_EXPERIMENTS_HOME)/software/art/art_illumina
	$$(QTIP_EXPERIMENTS_HOME)/software/art/art_illumina -sam -na -rs $6 -l $3 -f $5 -i $2 -o .$$@
	python $(QTIP_EXPERIMENTS_HOME)/bin/art_convert.py --in1 .$$(@).fq --sam .$$(@).sam --out1 .$$(@).final.fq
	rm -f .$$(@).fq .$$(@).sam
	head -n `python -c 'print(4 * $4)'` .$$(@).final.fq | gzip -c > $$@
	rm -f .$$(@).final.fq
	$$(QTIP_EXPERIMENTS_HOME)/software/art/art_illumina 2>&1 > $$@.version || true

endef

# Use art to generate paired-end Illumina-like reads.  Macro takes
# parameters: (1) batch name, (2) reference genome, (3) read length in
# batch, (4) target # pairs, (5) average fold-coverage, (6) mean fragment
# length, (7) stddev for fragment length, (8) pseudo-random seed.
define art_ill_pair_reads

# Generate paired-end reads
r1_art_$1.fq.gz: $(FA) $$(FA) $$(QTIP_EXPERIMENTS_HOME)/software/art/art_illumina
	$$(QTIP_EXPERIMENTS_HOME)/software/art/art_illumina -sam -na -p -rs $8 -m $6 -s $7 -l $3 -f $5 -i $2 -o .$$@
	python $(QTIP_EXPERIMENTS_HOME)/bin/art_convert.py --in1 .$$(@)1.fq --in2 .$$(@)2.fq --sam .$$(@).sam --out1 .$$(@).final1.fq --out2 .$$(@).final2.fq
	rm -f .$$(@)1.fq .$$(@)2.fq .$$(@).sam
	head -n `python -c 'print(4 * $4)'` .$$(@).final1.fq | gzip -c > $$@
	head -n `python -c 'print(4 * $4)'` .$$(@).final2.fq | gzip -c > $$(@:r1_%=r2_%)
	rm -f .$$(@).final1.fq .$$(@).final2.fq
	$$(QTIP_EXPERIMENTS_HOME)/software/art/art_illumina 2>&1 > $$@.version || true

endef

$(QTIP_EXPERIMENTS_HOME)/software/art/art_illumina:
	$(MAKE) -C $(QTIP_EXPERIMENTS_HOME)/software/art
