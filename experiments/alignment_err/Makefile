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

MASON=$(QSIM_EXPERIMENTS_HOME)/software/mason/mason

.PHONY: all
all: r0_mason_human_mixture_100_1M.fq \
     r0_mason_human_mixture_250_1M.fq \
     r1_mason_human_mixture_100_1M.fq \
     r1_mason_human_mixture_250_1M.fq

# 93% 1 target species (human or mouse)
#  6% 6 contaminants
#  1% 1 other mammal (either human or mouse, whichever isn't target)

r0_mason_human_mixture_100_1M.fq: reads
	gzip -dc ../simulated_reads/various_genomes/r0_mason_ill_hg_100_1M.fq.gz | head -n `expr 930000 \* 4` > $@
	gzip -dc ../simulated_reads/various_genomes/r0_mason_ill_mm_100_1M.fq.gz | head -n 40000 >> $@
	cat r0_mason_a_laidlawii_100_10k.fq >> $@
	cat r0_mason_m_fermentans_100_10k.fq >> $@
	cat r0_mason_m_hyorhinis_100_10k.fq >> $@
	cat r0_mason_m_hominis_100_10k.fq >> $@
	cat r0_mason_m_globosa_100_10k.fq >> $@
	cat r0_mason_p_acnes_100_10k.fq >> $@

r0_mason_human_mixture_250_1M.fq: reads
	gzip -dc ../simulated_reads/various_genomes/r0_mason_ill_hg_250_1M.fq.gz | head -n `expr 930000 \* 4` > $@
	gzip -dc ../simulated_reads/various_genomes/r0_mason_ill_mm_250_1M.fq.gz | head -n 40000 >> $@
	cat r0_mason_a_laidlawii_250_10k.fq >> $@
	cat r0_mason_m_fermentans_250_10k.fq >> $@
	cat r0_mason_m_hyorhinis_250_10k.fq >> $@
	cat r0_mason_m_hominis_250_10k.fq >> $@
	cat r0_mason_m_globosa_250_10k.fq >> $@
	cat r0_mason_p_acnes_250_10k.fq >> $@

r1_mason_human_mixture_100_1M.fq: reads
	gzip -dc ../simulated_reads/various_genomes/r1_mason_ill_hg_100_1M.fq.gz | head -n `expr 930000 \* 4` > $@
	gzip -dc ../simulated_reads/various_genomes/r2_mason_ill_hg_100_1M.fq.gz | head -n `expr 930000 \* 4` > $(@:r1_%=r2_%)
	gzip -dc ../simulated_reads/various_genomes/r1_mason_ill_mm_100_1M.fq.gz | head -n 40000 >> $@
	gzip -dc ../simulated_reads/various_genomes/r2_mason_ill_mm_100_1M.fq.gz | head -n 40000 >> $(@:r1_%=r2_%)
	cat r1_mason_a_laidlawii_100_10k.fq >> $@
	cat r1_mason_m_fermentans_100_10k.fq >> $@
	cat r1_mason_m_hyorhinis_100_10k.fq >> $@
	cat r1_mason_m_hominis_100_10k.fq >> $@
	cat r1_mason_m_globosa_100_10k.fq >> $@
	cat r1_mason_p_acnes_100_10k.fq >> $@
	cat r2_mason_a_laidlawii_100_10k.fq >> $(@:r1_%=r2_%)
	cat r2_mason_m_fermentans_100_10k.fq >> $(@:r1_%=r2_%)
	cat r2_mason_m_hyorhinis_100_10k.fq >> $(@:r1_%=r2_%)
	cat r2_mason_m_hominis_100_10k.fq >> $(@:r1_%=r2_%)
	cat r2_mason_m_globosa_100_10k.fq >> $(@:r1_%=r2_%)
	cat r2_mason_p_acnes_100_10k.fq >> $(@:r1_%=r2_%)

r1_mason_human_mixture_250_1M.fq: reads
	gzip -dc ../simulated_reads/various_genomes/r1_mason_ill_hg_250_1M.fq.gz | head -n `expr 930000 \* 4` > $@
	gzip -dc ../simulated_reads/various_genomes/r2_mason_ill_hg_250_1M.fq.gz | head -n `expr 930000 \* 4` > $(@:r1_%=r2_%)
	gzip -dc ../simulated_reads/various_genomes/r1_mason_ill_mm_250_1M.fq.gz | head -n 40000 >> $@
	gzip -dc ../simulated_reads/various_genomes/r2_mason_ill_mm_250_1M.fq.gz | head -n 40000 >> $(@:r1_%=r2_%)
	cat r1_mason_a_laidlawii_250_10k.fq >> $@
	cat r1_mason_m_fermentans_250_10k.fq >> $@
	cat r1_mason_m_hyorhinis_250_10k.fq >> $@
	cat r1_mason_m_hominis_250_10k.fq >> $@
	cat r1_mason_m_globosa_250_10k.fq >> $@
	cat r1_mason_p_acnes_250_10k.fq >> $@
	cat r2_mason_a_laidlawii_250_10k.fq >> $(@:r1_%=r2_%)
	cat r2_mason_m_fermentans_250_10k.fq >> $(@:r1_%=r2_%)
	cat r2_mason_m_hyorhinis_250_10k.fq >> $(@:r1_%=r2_%)
	cat r2_mason_m_hominis_250_10k.fq >> $(@:r1_%=r2_%)
	cat r2_mason_m_globosa_250_10k.fq >> $(@:r1_%=r2_%)
	cat r2_mason_p_acnes_250_10k.fq >> $(@:r1_%=r2_%)

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

define reads

r0_mason_$(1)_100_10k.fq: $(2) $(MASON)
	cp $$< $$(<:%.fna=%.fa)
	$$(MASON) illumina -hn 2 -i -s 453 -sq -n 100 -N 10000 -o .$$@.fq $$(<:%.fna=%.fa)
	rm -f .$$@.fq.sam
	python $$(QSIM_EXPERIMENTS_HOME)/bin/mason_convert.py --in1 .$$@.fq --out1 $$@
	rm -f .$$@.fq

r0_mason_$(1)_250_10k.fq: $(2) $(MASON)
	cp $$< $$(<:%.fna=%.fa)
	$$(MASON) illumina -hn 2 -i -s 278 -sq -n 250 -N 10000 -o .$$@.fq $$(<:%.fna=%.fa)
	rm -f .$$@.fq.sam
	python $$(QSIM_EXPERIMENTS_HOME)/bin/mason_convert.py --in1 .$$@.fq --out1 $$@
	rm -f .$$@.fq

r1_mason_$(1)_100_10k.fq: $(2) $(MASON)
	cp $$< $$(<:%.fna=%.fa)
	$$(MASON) illumina -hn 2 -i -s 645 -sq -mp -rn 2 -ll 250 -le 100 -n 100 -N 10000 -o .$$@.fq $$(<:%.fna=%.fa)
	rm -f .$$@.fq.sam
	python $$(QSIM_EXPERIMENTS_HOME)/bin/mason_convert.py --in1 .$$(@)_1.fq --in2 .$$(@)_2.fq --out1 .$$(@)_final_1.fq --out2 .$$(@)_final_2.fq
	rm -f .$$(@)_1.fq .$$(@)_2.fq
	mv .$$(@)_final_1.fq $$@
	mv .$$(@)_final_2.fq $$(@:r1_%=r2_%)

r1_mason_$(1)_250_10k.fq: $(2) $(MASON)
	cp $$< $$(<:%.fna=%.fa)
	$$(MASON) illumina -hn 2 -i -s 502 -sq -mp -rn 2 -ll 500 -le 200 -n 250 -N 10000 -o .$$@.fq $$(<:%.fna=%.fa)
	rm -f .$$@.fq.sam
	python $$(QSIM_EXPERIMENTS_HOME)/bin/mason_convert.py --in1 .$$(@)_1.fq --in2 .$$(@)_2.fq --out1 .$$(@)_final_1.fq --out2 .$$(@)_final_2.fq
	rm -f .$$(@)_1.fq .$$(@)_2.fq
	mv .$$(@)_final_1.fq $$@
	mv .$$(@)_final_2.fq $$(@:r1_%=r2_%)

endef

$(eval $(call reads,p_acnes,GCA_001481615.1_ASM148161v1_genomic.fna))
$(eval $(call reads,m_globosa,GCA_001264815.1_ASM126481v1_genomic.fna))
$(eval $(call reads,m_hominis,GCA_001063305.1_ASM106330v1_genomic.fna))
$(eval $(call reads,m_hyorhinis,GCA_000496815.1_ASM49681v1_genomic.fna))
$(eval $(call reads,m_fermentans,GCA_000209735.1_ASM20973v1_genomic.fna))
$(eval $(call reads,a_laidlawii,GCA_000018785.1_ASM1878v1_genomic.fna))

.PHONY: clean
clean:
	rm -f .*.fq *.fq *.fa
