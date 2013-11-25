# art_illumina [options] -i <DNA_reference_file> -l <read_length> -f <fold_coverage> -o <outFile_prefix>
#
# Relevant options:
#  -mp  --matepair indicate a mate-pair read simulation
#  -p   --paired   indicate a paired-end read simulation
#       NOTE: art will automatically switch to a mate-pair read simulation if the given mean fragment size >= 2000
#  -l   --len      the length of reads to be simulated
#  -f   --fcov     the fold of read coverage to be simulated
#  -m   --mflen    the mean size of DNA fragments for paired-end simulations
#  -s   --sdev     the standard deviation of DNA fragment size for paired-end simulations.
#  -na  --noALN    do not output ALN alignment file
#  -rs  --rndSend  the seed for random number generator (default: system time in second)

# Use wgsim to generate unpaired Illumina-like reads.  Macro takes
# parameters: (1) batch name, (2) reference genome, (3) read length in
# batch, (4) # reads in batch, (5) pseudo-random seed.
define art_ill_unp_reads

# Generate unpaired reads
r0_art_$1.fq.gz: $$(FA) $$(TS_HOME)/software/art/art_illumina
	$$(TS_HOME)/software/art/art_illumina -rs $5 -l $3 -f $4 -na -i $2 -o $$(@:%.fq.gz=%)
	gzip $$(@:%.fq.gz=%.fq)

endef

# Use wgsim to generate paired-end Illumina-like reads.  Macro takes
# parameters: (1) batch name, (2) reference genome, (3) read length in
# batch, (4) average fold-coverage, (5) mean fragment length, (6)
# stddev for fragment length, (7) pseudo-random seed.
define art_ill_pair_reads

# Generate paired-end reads
r1_art_$1.fq.gz: $(FA) $$(FA) $$(TS_HOME)/software/art/art_illumina
	$$(TS_HOME)/software/art/art_illumina -rs $7 -m $5 -s $6 -l $3 -f $4 -na -i $2 -o .tmp
	gzip -c .tmp1.fq > $$@
	rm -f .tmp1.fq
	gzip -c .tmp2.fq > $$(@:r1_%=r2_%)
	rm -f .tmp2.fq

endef

$(TS_HOME)/software/art/art_illumina:
	$(MAKE) -C $(TS_HOME)/software/art
