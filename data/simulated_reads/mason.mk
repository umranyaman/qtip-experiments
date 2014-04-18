# Relevant Mason parameters:
# -N   # reads
# -i   Include additional read information in reads file
# -sq  Simulate qualities
# -hn  # haplotypes to simulate
# -hs  haplotype SNP rate (0.001)
# -hi  haplotype indel rate (0.001)
# -hm  Haplotype indel size min (1)
# -hM  Haplotype indel size max (6)
# Illumina:
# -n   read length for illumina
# 454:
# -nu  Use uniform read length distribution (default: normal)
# -nm  Read length mean
# -ne  Read length error (stddev for normal, interval for uniform)

# Use Mason to generate unpaired Illumina-like reads.  Macro takes
# parameters: (1) batch name, (2) reference genome, (3) read length in
# batch, (4) # reads in batch, (5) pseudo-random seed.
define mason_ill_unp_reads

# Generate unpaired reads
r0_mason_$1.fq.gz: $$(FA) $$(TS_HOME)/software/mason/mason
	$$(TS_HOME)/software/mason/mason illumina -hn 2 -i -s $5 -sq -n $3 -N $4 -o .$$@.fq $2
	rm -f .$$@.fq.sam
	python $$(TS_HOME)/bin/mason_convert.py --in1 .$$@.fq --out1 $$(@:%.fq.gz=%.fq)
	rm -f .$$@.fq
	gzip $$(@:%.fq.gz=%.fq)
	$$(TS_HOME)/software/mason/mason illumina --version > $$@.version

endef

# Use Mason to generate paired-end Illumina-like reads.  Macro takes
# parameters: (1) batch name, (2) reference genome, (3) read length in
# batch, (4) # reads in batch, (5) mean fragment length, (6) stddev for
# fragment length, (7) pseudo-random seed
define mason_ill_pair_reads

# Generate paired-end reads
r1_mason_$1.fq.gz: $(FA) $$(TS_HOME)/software/mason/mason
	$$(TS_HOME)/software/mason/mason illumina -hn 2 -i -s $7 -sq -mp -rn 2 -ll $5 -le $6 -n $3 -N $4 -o .$$@.fq $2
	rm -f .$$@.fq.sam
	python $$(TS_HOME)/bin/mason_convert.py --in1 .$$(@)_1.fq --in2 .$$(@)_2.fq --out1 .$$(@)_final_1.fq --out2 .$$(@)_final_2.fq
	rm -f .$$(@)_1.fq .$$(@)_2.fq
	gzip -c .$$(@)_final_1.fq > $$@
	rm -f .$$(@)_final_1.fq
	gzip -c .$$(@)_final_2.fq > $$(@:r1_%=r2_%)
	rm -f .$$(@)_final_2.fq
	$$(TS_HOME)/software/mason/mason illumina --version > $$@.version

endef

# Use Mason to generate unpaired 454-like reads of various lengths.
# Macro takes parameters: (1) batch name, (2) reference genome, (3)
# read length mean, (4) read length error, (5) # reads in batch, (6)
# pseudo-random seed
define mason_fff_unp_reads

# Generate unpaired reads
r0_mason_$1.fq.gz: $$(FA) $$(TS_HOME)/software/mason/mason
	$$(TS_HOME)/software/mason/mason 454 -hn 2 -i -s $6 -sq -N $5 -nm $3 -ne $4 -nu -o $$(@:%.fq.gz=%.fq) $2
	rm -f $$(@:%.gz=%.sam)
	gzip $$(@:%.fq.gz=%.fq)
	$$(TS_HOME)/software/mason/mason 454 --version > $$@.version

endef

$(TS_HOME)/software/mason/mason:
	$(MAKE) -C $(TS_HOME)/software/mason
