# Program: wgsim (short read simulator)
# Version: 0.3.1-r13
# Contact: Heng Li <lh3@sanger.ac.uk>
#
# Usage:   wgsim [options] <in.ref.fa> <out.read1.fq> <out.read2.fq>
#
# Options: -e FLOAT      base error rate [0.020]
#          -d INT        outer distance between the two ends [500]
#          -s INT        standard deviation [50]
#          -N INT        number of read pairs [1000000]
#          -1 INT        length of the first read [70]
#          -2 INT        length of the second read [70]
#          -r FLOAT      rate of mutations [0.0010]
#          -R FLOAT      fraction of indels [0.15]
#          -X FLOAT      probability an indel is extended [0.30]
#          -S INT        seed for random generator [-1]
#          -A FLOAT      disgard if the fraction of ambiguous bases higher than FLOAT [0.05]
#          -h            haplotype mode

# Use wgsim to generate unpaired Illumina-like reads.  Macro takes
# parameters: (1) batch name, (2) reference genome, (3) read length in
# batch, (4) # reads in batch, (5) pseudo-random seed.
define wgsim_ill_unp_reads

# Generate unpaired reads
r0_wgsim_$1.fq.gz: $$(FA) $$(QSIM_EXPERIMENTS_HOME)/software/wgsim/wgsim
	$$(QSIM_EXPERIMENTS_HOME)/software/wgsim/wgsim -S $5 -1 $3 -2 $3 -N $4 $2 $$(@:%.fq.gz=%.fq) .tmp.$$@
	rm -f .tmp.$$@
	gzip $$(@:%.fq.gz=%.fq)
	$$(QSIM_EXPERIMENTS_HOME)/software/wgsim/wgsim > $$@.version 2>&1 || true

endef

# Use wgsim to generate paired-end Illumina-like reads.  Macro takes
# parameters: (1) batch name, (2) reference genome, (3) read length in
# batch, (4) # reads in batch, (5) mean fragment length, (6) stddev for
# fragment length, (7) pseudo-random seed.
define wgsim_ill_pair_reads

# Generate paired-end reads
r1_wgsim_$1.fq.gz: $(FA) $$(QSIM_EXPERIMENTS_HOME)/software/wgsim/wgsim
	$$(QSIM_EXPERIMENTS_HOME)/software/wgsim/wgsim -S $7 -d $5 -s $6 -1 $3 -2 $3 -N $4 $2 .tmp_1.$$@.fq .tmp_2.$$@.fq
	gzip -c .tmp_1.$$@.fq > $$@
	rm -f .tmp_1.$$@.fq
	gzip -c .tmp_2.$$@.fq > $$(@:r1_%=r2_%)
	rm -f .tmp_2.$$@.fq
	$$(QSIM_EXPERIMENTS_HOME)/software/wgsim/wgsim 2>&1 > $$@.version || true

endef

$(QSIM_EXPERIMENTS_HOME)/software/wgsim/wgsim:
	$(MAKE) -C $(QSIM_EXPERIMENTS_HOME)/software/wgsim
