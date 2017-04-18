# Relevant Mason parameters:
# -N   # reads
# -s   --seed INT random seed
# -i   Include additional read information in reads file
# -sq  Simulate qualities
# -hn  # haplotypes to simulate
# -hs  haplotype SNP rate (0.001)
# -hi  haplotype indel rate (0.001)
# -hm  Haplotype indel size min (1)
# -hM  Haplotype indel size max (6)
# Illumina:
# -n   read length for illumina
# -ll, --library-length-mean NUM
# -le, --library-length-error NUM
# 454:
# -nu  Use uniform read length distribution (default: normal)
# -nm  Read length mean
# -ne  Read length error (stddev for normal, interval for uniform)

ifndef QTIP_HOME
$(error QTIP_HOME is not set)
endif
ifndef QTIP_EXPERIMENTS_HOME
$(error QTIP_EXPERIMENTS_HOME is not set)
endif

MASON=$(QTIP_EXPERIMENTS_HOME)/software/mason/mason
BT2=$(QTIP_HOME)/software/bowtie2/bowtie2
NTHREADS=24
SHELL=/bin/bash
ASM_HG38=hg38
ASM_HG19=hg19
ASM_MOUSE=mm10
ASM_HG38_DECOY=$(QTIP_EXPERIMENTS_HOME)/experiments/refs/hg38_decoy.fa
ASM_HG19_DECOY=$(QTIP_EXPERIMENTS_HOME)/experiments/refs/hg19_decoy.fa
ASM_HG38_CHM1=$(QTIP_EXPERIMENTS_HOME)/experiments/refs/hg38_with_chm1
ASM_HG19_CHM1=$(QTIP_EXPERIMENTS_HOME)/experiments/refs/hg19_with_chm1

.PHONY: all
all: corstats

.PHONY: corstats
corstats: r0_mason_hg38_mixture_100.corstats \
          r0_mason_hg38_mixture_250.corstats \
          r1_mason_hg38_mixture_100.corstats \
          r1_mason_hg38_mixture_250.corstats \
          r0_mason_hg38chm1_mixture_100.corstats \
          r0_mason_hg38chm1_mixture_250.corstats \
          r1_mason_hg38chm1_mixture_100.corstats \
          r1_mason_hg38chm1_mixture_250.corstats \
          r0_mason_mouse_mixture_100.corstats \
          r0_mason_mouse_mixture_250.corstats \
          r1_mason_mouse_mixture_100.corstats \
          r1_mason_mouse_mixture_250.corstats

r0_mason_hg38_%.corstats: r0_mason_hg38_%.sam
	pypy evaluate.py human < $< > $@

r0_mason_hg38chm1_%.corstats: r0_mason_hg38_%.sam
	pypy evaluate.py human < $< > $@

r0_mason_mouse_%.corstats: r0_mason_mouse_%.sam
	pypy evaluate.py mouse < $< > $@

r1_mason_hg38_mixture_100.corstats: r1_mason_hg38_mixture_100.sam
	pypy evaluate.py human < $< > $@

r1_mason_hg38_mixture_250.corstats: r1_mason_hg38_mixture_250.sam
	pypy evaluate.py human < $< > $@

r1_mason_hg38chm1_mixture_100.corstats: r1_mason_hg38chm1_mixture_100.sam
	pypy evaluate.py human < $< > $@

r1_mason_hg38chm1_mixture_250.corstats: r1_mason_hg38chm1_mixture_250.sam
	pypy evaluate.py human < $< > $@

r1_mason_mouse_mixture_100.corstats: r1_mason_mouse_mixture_100.sam
	pypy evaluate.py mouse < $< > $@

r1_mason_mouse_mixture_250.corstats: r1_mason_mouse_mixture_250.sam
	pypy evaluate.py mouse < $< > $@

.PHONY: alignments
alignments: r0_mason_hg38_mixture_100.sam \
            r0_mason_hg38_mixture_250.sam \
            r1_mason_hg38_mixture_100.sam \
            r1_mason_hg38_mixture_250.sam \
            r0_mason_hg38chm1_mixture_100.sam \
            r0_mason_hg38chm1_mixture_250.sam \
            r1_mason_hg38chm1_mixture_100.sam \
            r1_mason_hg38chm1_mixture_250.sam \
            r0_mason_mouse_mixture_100.sam \
            r0_mason_mouse_mixture_250.sam \
            r1_mason_mouse_mixture_100.sam \
            r1_mason_mouse_mixture_250.sam

# 93% 1 target species (human or mouse)
#  1% 1 other mammal (either human or mouse, whichever isn't target)
#  6% 6 contaminants

# The target species will be the reference genome.  After aligning the mixture
# of reads to the reference, we want to compile
#
# 1. Type 1 errors: # reads from a non-target aligning to target
# 2. Type 2 errors: # reads simulated from target that failed to align
# 3. Type 3 errors: # reads that aligned to wrong place on target
#
# As such, the only read names that are important are the names of the reads
# simulated from the target.  That's the reason for the $(RENAME commands
# below.

RENAME=awk '{ln += 1; if(ln % 4 == 1) { n += 1; print "@r"n } else {print}}'

r0_mason_hg38_%.sam: r0_mason_hg38_%.fq
	$(BT2) \
	    -x $(QTIP_EXPERIMENTS_HOME)/experiments/refs/$(ASM_HG38).fa \
	    --sam-no-qname-trunc \
	    -p $(NTHREADS) \
	    -U $< \
	    -t \
	    -S $@ 2> >(tee $@.alstats)

r0_mason_hg38chm1_%.sam: r0_mason_hg38chm1_%.fq
	$(BT2) \
	    -x $(QTIP_EXPERIMENTS_HOME)/experiments/refs/$(ASM_HG38).fa \
	    --sam-no-qname-trunc \
	    -p $(NTHREADS) \
	    -U $< \
	    -t \
	    -S $@ 2> >(tee $@.alstats)

r0_mason_mouse_%.sam: r0_mason_mouse_%.fq
	$(BT2) \
	    -x $(QTIP_EXPERIMENTS_HOME)/experiments/refs/$(ASM_MOUSE).fa \
	    --sam-no-qname-trunc \
	    -p $(NTHREADS) \
	    -U $< \
	    -t \
	    -S $@ 2> >(tee $@.alstats)

r1_mason_hg38_mixture_100.sam: r1_mason_hg38_mixture_100.fq
	$(BT2) \
	    -x $(QTIP_EXPERIMENTS_HOME)/experiments/refs/$(ASM_HG38).fa \
	    --sam-no-qname-trunc \
	    -p $(NTHREADS) \
	    -1 $< \
	    -2 $(<:r1%=r2%) \
	    -I 200 -X 400 \
	    -t \
	    -S $@ 2> >(tee $@.alstats)

r1_mason_hg38_mixture_250.sam: r1_mason_hg38_mixture_250.fq
	$(BT2) \
	    -x $(QTIP_EXPERIMENTS_HOME)/experiments/refs/$(ASM_HG38).fa \
	    --sam-no-qname-trunc \
	    -p $(NTHREADS) \
	    -1 $< \
	    -2 $(<:r1%=r2%) \
	    -I 500 -X 1000 \
	    -t \
	    -S $@ 2> >(tee $@.alstats)

r1_mason_hg38chm1_mixture_100.sam: r1_mason_hg38chm1_mixture_100.fq
	$(BT2) \
	    -x $(QTIP_EXPERIMENTS_HOME)/experiments/refs/$(ASM_HG38).fa \
	    --sam-no-qname-trunc \
	    -p $(NTHREADS) \
	    -1 $< \
	    -2 $(<:r1%=r2%) \
	    -I 200 -X 400 \
	    -t \
	    -S $@ 2> >(tee $@.alstats)

r1_mason_hg38chm1_mixture_250.sam: r1_mason_hg38chm1_mixture_250.fq
	$(BT2) \
	    -x $(QTIP_EXPERIMENTS_HOME)/experiments/refs/$(ASM_HG38).fa \
	    --sam-no-qname-trunc \
	    -p $(NTHREADS) \
	    -1 $< \
	    -2 $(<:r1%=r2%) \
	    -I 500 -X 1000 \
	    -t \
	    -S $@ 2> >(tee $@.alstats)

r1_mason_mouse_mixture_100.sam: r1_mason_mouse_mixture_100.fq
	$(BT2) \
	    -x $(QTIP_EXPERIMENTS_HOME)/experiments/refs/$(ASM_MOUSE).fa \
	    --sam-no-qname-trunc \
	    -p $(NTHREADS) \
	    -1 $< \
	    -2 $(<:r1%=r2%) \
	    -I 200 -X 400 \
	    -t \
	    -S $@ 2> >(tee $@.alstats)

r1_mason_mouse_mixture_250.sam: r1_mason_mouse_mixture_250.fq
	$(BT2) \
	    -x $(QTIP_EXPERIMENTS_HOME)/experiments/refs/$(ASM_MOUSE).fa \
	    --sam-no-qname-trunc \
	    -p $(NTHREADS) \
	    -1 $< \
	    -2 $(<:r1%=r2%) \
	    -I 500 -X 1000 \
	    -t \
	    -S $@ 2> >(tee $@.alstats)

# What if I want these experiments to include CHM1 reads?  I could do this
# but point to the various_genomes_chm1 directory instead of various_genomes.
# Should I do both?

READS_UNP_100=../simulated_reads/various_genomes/r0_mason_ill_hg38_100.fq.gz \
              ../simulated_reads/various_genomes/r0_mason_ill_mm_100.fq.gz \
              r0_mason_a_laidlawii_100_10k.fq \
              r0_mason_m_fermentans_100_10k.fq \
              r0_mason_m_hyorhinis_100_10k.fq \
              r0_mason_m_hominis_100_10k.fq \
              r0_mason_m_globosa_100_10k.fq \
              r0_mason_p_acnes_100_10k.fq

READS_UNP_CHM1_100=../simulated_reads/various_genomes_chm1/r0_mason_ill_hg38_100.fq.gz \
              ../simulated_reads/various_genomes/r0_mason_ill_mm_100.fq.gz \
              r0_mason_a_laidlawii_100_10k.fq \
              r0_mason_m_fermentans_100_10k.fq \
              r0_mason_m_hyorhinis_100_10k.fq \
              r0_mason_m_hominis_100_10k.fq \
              r0_mason_m_globosa_100_10k.fq \
              r0_mason_p_acnes_100_10k.fq

READS_UNP_250=../simulated_reads/various_genomes/r0_mason_ill_hg38_250.fq.gz \
              ../simulated_reads/various_genomes/r0_mason_ill_mm_250.fq.gz \
              r0_mason_a_laidlawii_250_10k.fq \
              r0_mason_m_fermentans_250_10k.fq \
              r0_mason_m_hyorhinis_250_10k.fq \
              r0_mason_m_hominis_250_10k.fq \
              r0_mason_m_globosa_250_10k.fq \
              r0_mason_p_acnes_250_10k.fq

READS_UNP_CHM1_250=../simulated_reads/various_genomes_chm1/r0_mason_ill_hg38_250.fq.gz \
              ../simulated_reads/various_genomes/r0_mason_ill_mm_250.fq.gz \
              r0_mason_a_laidlawii_250_10k.fq \
              r0_mason_m_fermentans_250_10k.fq \
              r0_mason_m_hyorhinis_250_10k.fq \
              r0_mason_m_hominis_250_10k.fq \
              r0_mason_m_globosa_250_10k.fq \
              r0_mason_p_acnes_250_10k.fq

READS_PAIR_100=../simulated_reads/various_genomes/r1_mason_ill_hg38_100.fq.gz \
               ../simulated_reads/various_genomes/r1_mason_ill_mm_100.fq.gz \
               r1_mason_a_laidlawii_100_10k.fq \
               r1_mason_m_fermentans_100_10k.fq \
               r1_mason_m_hyorhinis_100_10k.fq \
               r1_mason_m_hominis_100_10k.fq \
               r1_mason_m_globosa_100_10k.fq \
               r1_mason_p_acnes_100_10k.fq \
               ../simulated_reads/various_genomes/r2_mason_ill_hg38_100.fq.gz \
               ../simulated_reads/various_genomes/r2_mason_ill_mm_100.fq.gz \
               r2_mason_a_laidlawii_100_10k.fq \
               r2_mason_m_fermentans_100_10k.fq \
               r2_mason_m_hyorhinis_100_10k.fq \
               r2_mason_m_hominis_100_10k.fq \
               r2_mason_m_globosa_100_10k.fq \
               r2_mason_p_acnes_100_10k.fq

READS_PAIR_CHM1_100=../simulated_reads/various_genomes_chm1/r1_mason_ill_hg38_100.fq.gz \
               ../simulated_reads/various_genomes/r1_mason_ill_mm_100.fq.gz \
               r1_mason_a_laidlawii_100_10k.fq \
               r1_mason_m_fermentans_100_10k.fq \
               r1_mason_m_hyorhinis_100_10k.fq \
               r1_mason_m_hominis_100_10k.fq \
               r1_mason_m_globosa_100_10k.fq \
               r1_mason_p_acnes_100_10k.fq \
               ../simulated_reads/various_genomes_chm1/r2_mason_ill_hg38_100.fq.gz \
               ../simulated_reads/various_genomes/r2_mason_ill_mm_100.fq.gz \
               r2_mason_a_laidlawii_100_10k.fq \
               r2_mason_m_fermentans_100_10k.fq \
               r2_mason_m_hyorhinis_100_10k.fq \
               r2_mason_m_hominis_100_10k.fq \
               r2_mason_m_globosa_100_10k.fq \
               r2_mason_p_acnes_100_10k.fq

READS_PAIR_250=../simulated_reads/various_genomes/r1_mason_ill_hg38_250.fq.gz \
               ../simulated_reads/various_genomes/r1_mason_ill_mm_250.fq.gz \
               r1_mason_a_laidlawii_250_10k.fq \
               r1_mason_m_fermentans_250_10k.fq \
               r1_mason_m_hyorhinis_250_10k.fq \
               r1_mason_m_hominis_250_10k.fq \
               r1_mason_m_globosa_250_10k.fq \
               r1_mason_p_acnes_250_10k.fq \
               ../simulated_reads/various_genomes/r2_mason_ill_hg38_250.fq.gz \
               ../simulated_reads/various_genomes/r2_mason_ill_mm_250.fq.gz \
               r2_mason_a_laidlawii_250_10k.fq \
               r2_mason_m_fermentans_250_10k.fq \
               r2_mason_m_hyorhinis_250_10k.fq \
               r2_mason_m_hominis_250_10k.fq \
               r2_mason_m_globosa_250_10k.fq \
               r2_mason_p_acnes_250_10k.fq

READS_PAIR_CHM1_250=../simulated_reads/various_genomes_chm1/r1_mason_ill_hg38_250.fq.gz \
               ../simulated_reads/various_genomes/r1_mason_ill_mm_250.fq.gz \
               r1_mason_a_laidlawii_250_10k.fq \
               r1_mason_m_fermentans_250_10k.fq \
               r1_mason_m_hyorhinis_250_10k.fq \
               r1_mason_m_hominis_250_10k.fq \
               r1_mason_m_globosa_250_10k.fq \
               r1_mason_p_acnes_250_10k.fq \
               ../simulated_reads/various_genomes_chm1/r2_mason_ill_hg38_250.fq.gz \
               ../simulated_reads/various_genomes/r2_mason_ill_mm_250.fq.gz \
               r2_mason_a_laidlawii_250_10k.fq \
               r2_mason_m_fermentans_250_10k.fq \
               r2_mason_m_hyorhinis_250_10k.fq \
               r2_mason_m_hominis_250_10k.fq \
               r2_mason_m_globosa_250_10k.fq \
               r2_mason_p_acnes_250_10k.fq

r0_mason_hg38_mixture_100.fq: $(READS_UNP_100)
	gzip -dc ../simulated_reads/various_genomes/r0_mason_ill_hg38_100.fq.gz | head -n `expr 930000 \* 4` > $@
	gzip -dc ../simulated_reads/various_genomes/r0_mason_ill_mm_100.fq.gz | head -n 40000 | $(RENAME) >> $@
	cat r0_mason_a_laidlawii_100_10k.fq | $(RENAME) >> $@
	cat r0_mason_m_fermentans_100_10k.fq | $(RENAME) >> $@
	cat r0_mason_m_hyorhinis_100_10k.fq | $(RENAME) >> $@
	cat r0_mason_m_hominis_100_10k.fq | $(RENAME) >> $@
	cat r0_mason_m_globosa_100_10k.fq | $(RENAME) >> $@
	cat r0_mason_p_acnes_100_10k.fq | $(RENAME) >> $@

r0_mason_hg38chm1_mixture_100.fq: $(READS_UNP_CHM1_100)
	gzip -dc ../simulated_reads/various_genomes_chm1/r0_mason_ill_hg38_100.fq.gz | head -n `expr 930000 \* 4` > $@
	gzip -dc ../simulated_reads/various_genomes/r0_mason_ill_mm_100.fq.gz | head -n 40000 | $(RENAME) >> $@
	cat r0_mason_a_laidlawii_100_10k.fq | $(RENAME) >> $@
	cat r0_mason_m_fermentans_100_10k.fq | $(RENAME) >> $@
	cat r0_mason_m_hyorhinis_100_10k.fq | $(RENAME) >> $@
	cat r0_mason_m_hominis_100_10k.fq | $(RENAME) >> $@
	cat r0_mason_m_globosa_100_10k.fq | $(RENAME) >> $@
	cat r0_mason_p_acnes_100_10k.fq | $(RENAME) >> $@

r0_mason_mouse_mixture_100.fq: $(READS_UNP_100)
	gzip -dc ../simulated_reads/various_genomes/r0_mason_ill_mm_100.fq.gz | head -n `expr 930000 \* 4` > $@
	gzip -dc ../simulated_reads/various_genomes/r0_mason_ill_hg38_100.fq.gz | head -n 40000 | $(RENAME) >> $@
	cat r0_mason_a_laidlawii_100_10k.fq | $(RENAME) >> $@
	cat r0_mason_m_fermentans_100_10k.fq | $(RENAME) >> $@
	cat r0_mason_m_hyorhinis_100_10k.fq | $(RENAME) >> $@
	cat r0_mason_m_hominis_100_10k.fq | $(RENAME) >> $@
	cat r0_mason_m_globosa_100_10k.fq | $(RENAME) >> $@
	cat r0_mason_p_acnes_100_10k.fq | $(RENAME) >> $@

r0_mason_hg38_mixture_250.fq: $(READS_UNP_250)
	gzip -dc ../simulated_reads/various_genomes/r0_mason_ill_hg38_250.fq.gz | head -n `expr 930000 \* 4` > $@
	gzip -dc ../simulated_reads/various_genomes/r0_mason_ill_mm_250.fq.gz | head -n 40000 | $(RENAME) >> $@
	cat r0_mason_a_laidlawii_250_10k.fq | $(RENAME) >> $@
	cat r0_mason_m_fermentans_250_10k.fq | $(RENAME) >> $@
	cat r0_mason_m_hyorhinis_250_10k.fq | $(RENAME) >> $@
	cat r0_mason_m_hominis_250_10k.fq | $(RENAME) >> $@
	cat r0_mason_m_globosa_250_10k.fq | $(RENAME) >> $@
	cat r0_mason_p_acnes_250_10k.fq | $(RENAME) >> $@

r0_mason_hg38chm1_mixture_250.fq: $(READS_UNP_CHM1_250)
	gzip -dc ../simulated_reads/various_genomes_chm1/r0_mason_ill_hg38_250.fq.gz | head -n `expr 930000 \* 4` > $@
	gzip -dc ../simulated_reads/various_genomes/r0_mason_ill_mm_250.fq.gz | head -n 40000 | $(RENAME) >> $@
	cat r0_mason_a_laidlawii_250_10k.fq | $(RENAME) >> $@
	cat r0_mason_m_fermentans_250_10k.fq | $(RENAME) >> $@
	cat r0_mason_m_hyorhinis_250_10k.fq | $(RENAME) >> $@
	cat r0_mason_m_hominis_250_10k.fq | $(RENAME) >> $@
	cat r0_mason_m_globosa_250_10k.fq | $(RENAME) >> $@
	cat r0_mason_p_acnes_250_10k.fq | $(RENAME) >> $@

r0_mason_mouse_mixture_250.fq: $(READS_UNP_250)
	gzip -dc ../simulated_reads/various_genomes/r0_mason_ill_mm_250.fq.gz | head -n `expr 930000 \* 4` > $@
	gzip -dc ../simulated_reads/various_genomes/r0_mason_ill_hg38_250.fq.gz | head -n 40000 | $(RENAME) >> $@
	cat r0_mason_a_laidlawii_250_10k.fq | $(RENAME) >> $@
	cat r0_mason_m_fermentans_250_10k.fq | $(RENAME) >> $@
	cat r0_mason_m_hyorhinis_250_10k.fq | $(RENAME) >> $@
	cat r0_mason_m_hominis_250_10k.fq | $(RENAME) >> $@
	cat r0_mason_m_globosa_250_10k.fq | $(RENAME) >> $@
	cat r0_mason_p_acnes_250_10k.fq | $(RENAME) >> $@

r1_mason_hg38_mixture_100.fq: $(READS_PAIR_100)
	gzip -dc ../simulated_reads/various_genomes/r1_mason_ill_hg38_100.fq.gz | head -n `expr 930000 \* 4` > $@
	gzip -dc ../simulated_reads/various_genomes/r2_mason_ill_hg38_100.fq.gz | head -n `expr 930000 \* 4` > $(@:r1_%=r2_%)
	gzip -dc ../simulated_reads/various_genomes/r1_mason_ill_mm_100.fq.gz | head -n 40000 | $(RENAME) >> $@
	gzip -dc ../simulated_reads/various_genomes/r2_mason_ill_mm_100.fq.gz | head -n 40000 | $(RENAME) >> $(@:r1_%=r2_%)
	cat r1_mason_a_laidlawii_100_10k.fq | $(RENAME) >> $@
	cat r1_mason_m_fermentans_100_10k.fq | $(RENAME) >> $@
	cat r1_mason_m_hyorhinis_100_10k.fq | $(RENAME) >> $@
	cat r1_mason_m_hominis_100_10k.fq | $(RENAME) >> $@
	cat r1_mason_m_globosa_100_10k.fq | $(RENAME) >> $@
	cat r1_mason_p_acnes_100_10k.fq | $(RENAME) >> $@
	cat r2_mason_a_laidlawii_100_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_m_fermentans_100_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_m_hyorhinis_100_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_m_hominis_100_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_m_globosa_100_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_p_acnes_100_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)

r1_mason_hg38chm1_mixture_100.fq: $(READS_PAIR_CHM1_100)
	gzip -dc ../simulated_reads/various_genomes_chm1/r1_mason_ill_hg38_100.fq.gz | head -n `expr 930000 \* 4` > $@
	gzip -dc ../simulated_reads/various_genomes_chm1/r2_mason_ill_hg38_100.fq.gz | head -n `expr 930000 \* 4` > $(@:r1_%=r2_%)
	gzip -dc ../simulated_reads/various_genomes/r1_mason_ill_mm_100.fq.gz | head -n 40000 | $(RENAME) >> $@
	gzip -dc ../simulated_reads/various_genomes/r2_mason_ill_mm_100.fq.gz | head -n 40000 | $(RENAME) >> $(@:r1_%=r2_%)
	cat r1_mason_a_laidlawii_100_10k.fq | $(RENAME) >> $@
	cat r1_mason_m_fermentans_100_10k.fq | $(RENAME) >> $@
	cat r1_mason_m_hyorhinis_100_10k.fq | $(RENAME) >> $@
	cat r1_mason_m_hominis_100_10k.fq | $(RENAME) >> $@
	cat r1_mason_m_globosa_100_10k.fq | $(RENAME) >> $@
	cat r1_mason_p_acnes_100_10k.fq | $(RENAME) >> $@
	cat r2_mason_a_laidlawii_100_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_m_fermentans_100_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_m_hyorhinis_100_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_m_hominis_100_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_m_globosa_100_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_p_acnes_100_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)

r1_mason_mouse_mixture_100.fq: $(READS_PAIR_100)
	gzip -dc ../simulated_reads/various_genomes/r1_mason_ill_mm_100.fq.gz | head -n `expr 930000 \* 4` > $@
	gzip -dc ../simulated_reads/various_genomes/r2_mason_ill_mm_100.fq.gz | head -n `expr 930000 \* 4` > $(@:r1_%=r2_%)
	gzip -dc ../simulated_reads/various_genomes/r1_mason_ill_hg38_100.fq.gz | head -n 40000 | $(RENAME) >> $@
	gzip -dc ../simulated_reads/various_genomes/r2_mason_ill_hg38_100.fq.gz | head -n 40000 | $(RENAME) >> $(@:r1_%=r2_%)
	cat r1_mason_a_laidlawii_100_10k.fq | $(RENAME) >> $@
	cat r1_mason_m_fermentans_100_10k.fq | $(RENAME) >> $@
	cat r1_mason_m_hyorhinis_100_10k.fq | $(RENAME) >> $@
	cat r1_mason_m_hominis_100_10k.fq | $(RENAME) >> $@
	cat r1_mason_m_globosa_100_10k.fq | $(RENAME) >> $@
	cat r1_mason_p_acnes_100_10k.fq | $(RENAME) >> $@
	cat r2_mason_a_laidlawii_100_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_m_fermentans_100_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_m_hyorhinis_100_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_m_hominis_100_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_m_globosa_100_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_p_acnes_100_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)

r1_mason_hg38_mixture_250.fq: $(READS_PAIR_250)
	gzip -dc ../simulated_reads/various_genomes/r1_mason_ill_hg38_250.fq.gz | head -n `expr 930000 \* 4` > $@
	gzip -dc ../simulated_reads/various_genomes/r2_mason_ill_hg38_250.fq.gz | head -n `expr 930000 \* 4` > $(@:r1_%=r2_%)
	gzip -dc ../simulated_reads/various_genomes/r1_mason_ill_mm_250.fq.gz | head -n 40000 | $(RENAME) >> $@
	gzip -dc ../simulated_reads/various_genomes/r2_mason_ill_mm_250.fq.gz | head -n 40000 | $(RENAME) >> $(@:r1_%=r2_%)
	cat r1_mason_a_laidlawii_250_10k.fq | $(RENAME) >> $@
	cat r1_mason_m_fermentans_250_10k.fq | $(RENAME) >> $@
	cat r1_mason_m_hyorhinis_250_10k.fq | $(RENAME) >> $@
	cat r1_mason_m_hominis_250_10k.fq | $(RENAME) >> $@
	cat r1_mason_m_globosa_250_10k.fq | $(RENAME) >> $@
	cat r1_mason_p_acnes_250_10k.fq | $(RENAME) >> $@
	cat r2_mason_a_laidlawii_250_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_m_fermentans_250_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_m_hyorhinis_250_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_m_hominis_250_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_m_globosa_250_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_p_acnes_250_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)

r1_mason_hg38chm1_mixture_250.fq: $(READS_PAIR_CHM1_250)
	gzip -dc ../simulated_reads/various_genomes_chm1/r1_mason_ill_hg38_250.fq.gz | head -n `expr 930000 \* 4` > $@
	gzip -dc ../simulated_reads/various_genomes_chm1/r2_mason_ill_hg38_250.fq.gz | head -n `expr 930000 \* 4` > $(@:r1_%=r2_%)
	gzip -dc ../simulated_reads/various_genomes/r1_mason_ill_mm_250.fq.gz | head -n 40000 | $(RENAME) >> $@
	gzip -dc ../simulated_reads/various_genomes/r2_mason_ill_mm_250.fq.gz | head -n 40000 | $(RENAME) >> $(@:r1_%=r2_%)
	cat r1_mason_a_laidlawii_250_10k.fq | $(RENAME) >> $@
	cat r1_mason_m_fermentans_250_10k.fq | $(RENAME) >> $@
	cat r1_mason_m_hyorhinis_250_10k.fq | $(RENAME) >> $@
	cat r1_mason_m_hominis_250_10k.fq | $(RENAME) >> $@
	cat r1_mason_m_globosa_250_10k.fq | $(RENAME) >> $@
	cat r1_mason_p_acnes_250_10k.fq | $(RENAME) >> $@
	cat r2_mason_a_laidlawii_250_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_m_fermentans_250_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_m_hyorhinis_250_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_m_hominis_250_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_m_globosa_250_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_p_acnes_250_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)

r1_mason_mouse_mixture_250.fq: $(READS_PAIR_250)
	gzip -dc ../simulated_reads/various_genomes/r1_mason_ill_mm_250.fq.gz | head -n `expr 930000 \* 4` > $@
	gzip -dc ../simulated_reads/various_genomes/r2_mason_ill_mm_250.fq.gz | head -n `expr 930000 \* 4` > $(@:r1_%=r2_%)
	gzip -dc ../simulated_reads/various_genomes/r1_mason_ill_hg38_250.fq.gz | head -n 40000 | $(RENAME) >> $@
	gzip -dc ../simulated_reads/various_genomes/r2_mason_ill_hg38_250.fq.gz | head -n 40000 | $(RENAME) >> $(@:r1_%=r2_%)
	cat r1_mason_a_laidlawii_250_10k.fq | $(RENAME) >> $@
	cat r1_mason_m_fermentans_250_10k.fq | $(RENAME) >> $@
	cat r1_mason_m_hyorhinis_250_10k.fq | $(RENAME) >> $@
	cat r1_mason_m_hominis_250_10k.fq | $(RENAME) >> $@
	cat r1_mason_m_globosa_250_10k.fq | $(RENAME) >> $@
	cat r1_mason_p_acnes_250_10k.fq | $(RENAME) >> $@
	cat r2_mason_a_laidlawii_250_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_m_fermentans_250_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_m_hyorhinis_250_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_m_hominis_250_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_m_globosa_250_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)
	cat r2_mason_p_acnes_250_10k.fq | $(RENAME) >> $(@:r1_%=r2_%)

.PHONY: reads
reads: r0_mason_a_laidlawii_100_10k.fq \
       r0_mason_a_laidlawii_250_10k.fq \
       r1_mason_a_laidlawii_100_10k.fq \
       r1_mason_a_laidlawii_250_10k.fq \
       r0_mason_m_fermentans_100_10k.fq \
       r0_mason_m_fermentans_250_10k.fq \
       r1_mason_m_fermentans_100_10k.fq \
       r1_mason_m_fermentans_250_10k.fq \
       r0_mason_m_hyorhinis_100_10k.fq \
       r0_mason_m_hyorhinis_250_10k.fq \
       r1_mason_m_hyorhinis_100_10k.fq \
       r1_mason_m_hyorhinis_250_10k.fq \
       r0_mason_m_hominis_100_10k.fq \
       r0_mason_m_hominis_250_10k.fq \
       r1_mason_m_hominis_100_10k.fq \
       r1_mason_m_hominis_250_10k.fq \
       r0_mason_m_globosa_100_10k.fq \
       r0_mason_m_globosa_250_10k.fq \
       r1_mason_m_globosa_100_10k.fq \
       r1_mason_m_globosa_250_10k.fq \
       r0_mason_p_acnes_100_10k.fq \
       r0_mason_p_acnes_250_10k.fq \
       r1_mason_p_acnes_100_10k.fq \
       r1_mason_p_acnes_250_10k.fq

.PHONY: decoy_reads
decoy_reads: \
  r0_mason_hg38_decoy_100.fq \
  r0_mason_hg38_decoy_250.fq \
  r0_mason_hg19_decoy_100.fq \
  r0_mason_hg19_decoy_250.fq \
  r1_mason_hg38_decoy_100.fq \
  r1_mason_hg38_decoy_250.fq \
  r1_mason_hg19_decoy_100.fq \
  r1_mason_hg19_decoy_250.fq

.PHONY: chm1_reads
chm1_reads: \
  r0_mason_hg38_chm1_100.fq \
  r0_mason_hg38_chm1_250.fq \
  r0_mason_hg19_chm1_100.fq \
  r0_mason_hg19_chm1_250.fq \
  r1_mason_hg38_chm1_100.fq \
  r1_mason_hg38_chm1_250.fq \
  r1_mason_hg19_chm1_100.fq \
  r1_mason_hg19_chm1_250.fq

# GRCh38 100 nt alignments (all 3)
r0_mason_hg38_%_100.sam: r0_mason_hg38_%_100.fq
	$(BT2) \
	    -x $(QTIP_EXPERIMENTS_HOME)/experiments/refs/$(ASM_HG38).fa \
	    --sam-no-qname-trunc \
	    -U $< \
	    -t \
	    -S $@ 2> >(tee $@.alstats)

# GRCh38-CHM1 100 nt reads
r0_mason_hg38_chm1_100.fq: $(ASM_HG38_CHM1) $(MASON)
	$(MASON) illumina -i -s 453 -sq -n 100 -N 10000 -o $@ $< 2> .$@.err | tee .$@.out

# GRCh38-decoy 100 nt reads
r0_mason_hg38_decoy_100.fq: $(ASM_HG38_DECOY) $(MASON)
	$(MASON) illumina -i -s 453 -sq -n 100 -N 10000 -o $@ $< 2> .$@.err | tee .$@.out

# GRCh38 250 nt alignments (all 3)
r0_mason_hg38_%_250.sam: r0_mason_hg38_%_250.fq
	$(BT2) \
	    -x $(QTIP_EXPERIMENTS_HOME)/experiments/refs/$(ASM_HG38).fa \
	    --sam-no-qname-trunc \
	    -U $< \
	    -t \
	    -S $@ 2> >(tee $@.alstats)

# GRCh38-CHM1 250 nt reads
r0_mason_hg38_chm1_250.fq: $(ASM_HG38_CHM1) $(MASON)
	$(MASON) illumina -i -s 278 -sq -n 250 -N 10000 -o $@ $< 2> .$@.err | tee .$@.out

# GRCh38-decoy 250 nt reads
r0_mason_hg38_decoy_250.fq: $(ASM_HG38_DECOY) $(MASON)
	$(MASON) illumina -i -s 278 -sq -n 250 -N 10000 -o $@ $< 2> .$@.err | tee .$@.out

# GRCh37 100 nt alignments (all 3)
r0_mason_hg19_%_250.sam: r0_mason_hg19_%_250.fq
	$(BT2) \
	    -x $(QTIP_EXPERIMENTS_HOME)/experiments/refs/$(ASM_HG19).fa \
	    --sam-no-qname-trunc \
	    -U $< \
	    -t \
	    -S $@ 2> >(tee $@.alstats)

# GRCh37-CHM1 100 nt reads
r0_mason_hg19_chm1_100.fq: $(ASM_HG19_CHM1) $(MASON)
	$(MASON) illumina -i -s 453 -sq -n 100 -N 10000 -o $@ $< 2> .$@.err | tee .$@.out

# GRCh37-decoy 100 nt reads
r0_mason_hg19_decoy_100.fq: $(ASM_HG19_DECOY) $(MASON)
	$(MASON) illumina -i -s 453 -sq -n 100 -N 10000 -o $@ $< 2> .$@.err | tee .$@.out

# GRCh37 250 nt alignments (all 3)
r0_mason_hg19_%_250.sam: r0_mason_hg19_%_250.fq
	$(BT2) \
	    -x $(QTIP_EXPERIMENTS_HOME)/experiments/refs/$(ASM_HG19).fa \
	    --sam-no-qname-trunc \
	    -U $< \
	    -t \
	    -S $@ 2> >(tee $@.alstats)

# GRCh37-CHM1 250 nt reads
r0_mason_hg19_chm1_250.fq: $(ASM_HG19_CHM1) $(MASON)
	$(MASON) illumina -i -s 278 -sq -n 250 -N 10000 -o $@ $< 2> .$@.err | tee .$@.out

# GRCh37-decoy 250 nt reads
r0_mason_hg19_decoy_250.fq: $(ASM_HG19_DECOY) $(MASON)
	$(MASON) illumina -i -s 278 -sq -n 250 -N 10000 -o $@ $< 2> .$@.err | tee .$@.out

# GRCh37 100 nt paired-end alignments (all 3)
r1_mason_hg38_%_100.sam: r1_mason_hg38_%_100.fq
	$(BT2) \
	    -x $(QTIP_EXPERIMENTS_HOME)/experiments/refs/$(ASM_HG38).fa \
	    --sam-no-qname-trunc \
	    -1 $< \
	    -2 $(<:r1%=r2%) \
	    -I 200 -X 400 \
	    -t \
	    -S $@ 2> >(tee $@.alstats)

# GRCh38-CHM1 100 nt paired-end reads
r1_mason_hg38_chm1_100.fq: $(ASM_HG38_CHM1) $(MASON)
	$(MASON) illumina -i -s 645 -sq -mp -rn 2 -ll 300 -le 100 -n 100 -N 10000 -o .$(@).fq $< 2> .$@.err | tee .$@.out
	mv .$(@)_1.fq $@
	mv .$(@)_2.fq $(@:r1_%=r2_%)

# GRCh38-decoy 100 nt paired-end reads
r1_mason_hg38_decoy_100.fq: $(ASM_HG38_DECOY) $(MASON)
	$(MASON) illumina -i -s 645 -sq -mp -rn 2 -ll 300 -le 100 -n 100 -N 10000 -o .$(@).fq $< 2> .$@.err | tee .$@.out
	mv .$(@)_1.fq $@
	mv .$(@)_2.fq $(@:r1_%=r2_%)

# GRCh38 250 nt paired-end alignments (all 3)
r1_mason_hg38_%_250.sam: r1_mason_hg38_%_250.fq
	$(BT2) \
	    -x $(QTIP_EXPERIMENTS_HOME)/experiments/refs/$(ASM_HG38).fa \
	    --sam-no-qname-trunc \
	    -1 $< \
	    -2 $(<:r1%=r2%) \
	    -I 500 -X 1000 \
	    -t \
	    -S $@ 2> >(tee $@.alstats)

# GRCh38-CHM1 250 nt paired-end reads
r1_mason_hg38_chm1_250.fq: $(ASM_HG38_CHM1) $(MASON)
	$(MASON) illumina -i -s 502 -sq -mp -rn 2 -ll 750 -le 250 -n 250 -N 10000 -o .$(@).fq $< 2> .$@.err | tee .$@.out
	mv .$(@)_1.fq $@
	mv .$(@)_2.fq $(@:r1_%=r2_%)

# GRCh38-decoy 250 nt paired-end reads
r1_mason_hg38_decoy_250.fq: $(ASM_HG38_DECOY) $(MASON)
	$(MASON) illumina -i -s 502 -sq -mp -rn 2 -ll 750 -le 250 -n 250 -N 10000 -o .$(@).fq $< 2> .$@.err | tee .$@.out
	mv .$(@)_1.fq $@
	mv .$(@)_2.fq $(@:r1_%=r2_%)

# GRCh37 100 nt paired-end alignments (all 3)
r1_mason_hg19_%_100.sam: r1_mason_hg19_%_100.fq
	$(BT2) \
	    -x $(QTIP_EXPERIMENTS_HOME)/experiments/refs/$(ASM_HG19).fa \
	    --sam-no-qname-trunc \
	    -1 $< \
	    -2 $(<:r1%=r2%) \
	    -I 200 -X 400 \
	    -t \
	    -S $@ 2> >(tee $@.alstats)

# GRCh37-CHM1 100 nt paired-end reads
r1_mason_hg19_chm1_100.fq: $(ASM_HG19_CHM1) $(MASON)
	$(MASON) illumina -i -s 645 -sq -mp -rn 2 -ll 300 -le 100 -n 100 -N 10000 -o .$(@).fq $< 2> .$@.err | tee .$@.out
	mv .$(@)_1.fq $@
	mv .$(@)_2.fq $(@:r1_%=r2_%)

# GRCh37-decoy 100 nt paired-end reads
r1_mason_hg19_decoy_100.fq: $(ASM_HG19_DECOY) $(MASON)
	$(MASON) illumina -i -s 645 -sq -mp -rn 2 -ll 300 -le 100 -n 100 -N 10000 -o .$(@).fq $< 2> .$@.err | tee .$@.out
	mv .$(@)_1.fq $@
	mv .$(@)_2.fq $(@:r1_%=r2_%)

# GRCh37 250 nt paired-end alignments (all 3)
r1_mason_hg19_%_250.sam: r1_mason_hg19_%_250.fq
	$(BT2) \
	    -x $(QTIP_EXPERIMENTS_HOME)/experiments/refs/$(ASM_HG19).fa \
	    --sam-no-qname-trunc \
	    -1 $< \
	    -2 $(<:r1%=r2%) \
	    -I 500 -X 1000 \
	    -t \
	    -S $@ 2> >(tee $@.alstats)

# GRCh37-CHM1 100 nt paired-end reads
r1_mason_hg19_chm1_250.fq: $(ASM_HG19_CHM1) $(MASON)
	$(MASON) illumina -i -s 502 -sq -mp -rn 2 -ll 750 -le 250 -n 250 -N 10000 -o .$(@).fq $< 2> .$@.err | tee .$@.out
	mv .$(@)_1.fq $@
	mv .$(@)_2.fq $(@:r1_%=r2_%)

# GRCh37-decoy 100 nt paired-end reads
r1_mason_hg19_decoy_250.fq: $(ASM_HG19_DECOY) $(MASON)
	$(MASON) illumina -i -s 502 -sq -mp -rn 2 -ll 750 -le 250 -n 250 -N 10000 -o .$(@).fq $< 2> .$@.err | tee .$@.out
	mv .$(@)_1.fq $@
	mv .$(@)_2.fq $(@:r1_%=r2_%)

define reads

r0_mason_$(1)_100_10k.fq: $(2) $(MASON)
	cp $$< $$(<:%.fna=%.fa)
	$$(MASON) illumina -hn 2 -i -s 453 -sq -n 100 -N 10000 -o .$$@.fq $$(<:%.fna=%.fa)
	mv .$$@.fq $$@

r0_mason_$(1)_250_10k.fq: $(2) $(MASON)
	cp $$< $$(<:%.fna=%.fa)
	$$(MASON) illumina -hn 2 -i -s 278 -sq -n 250 -N 10000 -o .$$@.fq $$(<:%.fna=%.fa)
	mv .$$@.fq $$@

r1_mason_$(1)_100_10k.fq: $(2) $(MASON)
	cp $$< $$(<:%.fna=%.fa)
	$$(MASON) illumina -hn 2 -i -s 645 -sq -mp -rn 2 -ll 300 -le 100 -n 100 -N 10000 -o .$$@.fq $$(<:%.fna=%.fa)
	mv .$$(@)_1.fq $$@
	mv .$$(@)_2.fq $$(@:r1_%=r2_%)

r1_mason_$(1)_250_10k.fq: $(2) $(MASON)
	cp $$< $$(<:%.fna=%.fa)
	$$(MASON) illumina -hn 2 -i -s 502 -sq -mp -rn 2 -ll 750 -le 250 -n 250 -N 10000 -o .$$@.fq $$(<:%.fna=%.fa)
	mv .$$(@)_1.fq $$@
	mv .$$(@)_2.fq $$(@:r1_%=r2_%)

endef

$(eval $(call reads,p_acnes,GCA_001481615.1_ASM148161v1_genomic.fna))
$(eval $(call reads,m_globosa,GCA_001264815.1_ASM126481v1_genomic.fna))
$(eval $(call reads,m_hominis,GCA_001063305.1_ASM106330v1_genomic.fna))
$(eval $(call reads,m_hyorhinis,GCA_000496815.1_ASM49681v1_genomic.fna))
$(eval $(call reads,m_fermentans,GCA_000209735.1_ASM20973v1_genomic.fna))
$(eval $(call reads,a_laidlawii,GCA_000018785.1_ASM1878v1_genomic.fna))

.PHONY: clean
clean:
	rm -f .*.fq *.fa *.fq *.fna.gz *.fna
