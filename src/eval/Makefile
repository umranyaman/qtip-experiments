TS_BWA_HOME=$(TS_HOME)/software/bwa
TS_SNAP_HOME=$(TS_HOME)/software/snap
TS_BT2_HOME=$(TS_HOME)/software/bowtie2

BT2_ARGS=--aligner bowtie2 --bt2-exe $(TS_BT2_HOME)/bowtie2-git/bowtie2 --index $(TS_INDEXES)/hg19.fa
BWA_ARGS=--aligner bwa-mem --bwa-exe $(TS_BWA_HOME)/bwa --index $(TS_INDEXES)/hg19.fa
SNAP_ARGS=--aligner snap --snap-exe $(TS_SNAP_HOME)/snap/snap-aligner --index $(TS_INDEXES)/hg19.fa.snap

ARGS=--ref $(TS_REFS)/hg19.fa --write-orig-mapq --write-precise-mapq --seed 346
ARGS_SAMP=--keep-intermediates --predict-for-training

.PHONY: all
all: ERR050082_1.excerpt.csv ERR050083_1.excerpt.csv

%.csv: %.bwa.unp.sam %.bt2.unp.sam %.snap.unp.sam %.bt2vs.unp.sam %.bt2vsl.unp.sam
	python multi_aligner.py \
	    --sam $^ \
	    --name bwamem bt2 snap bt2vs bt2vsl \
	    --tier 1 1 2 2 2 > $@

even_tiered_paired.csv: bwamem_paired.sam bt2_paired.sam snap_paired.sam
	python multi_aligner.py --sam $^ --tier 1 1 1 > $@

# #######
# BWA-MEM
# #######

# BWA-MEM paired-end overall
%.bwa.pair.sam: %.fastq
	python ../ts.py \
		$(BWA_ARGS) $(ARGS) \
		--m1 $< \
		--m2 $(<:_1.excerpt.fastq=_2.excerpt.fastq) \
		--vanilla-output $@ \
		-- -t 8

# BWA-MEM paired-end samples
%.bwa.pair.sam.0: %.fastq.0
	for i in `seq 0 9`; do \
		python ../ts.py \
			$(BWA_ARGS) $(ARGS) $(ARGS_SAMP) \
			--m1 $(<:.0=.$$i) \
			--m2 $(<:_1.excerpt.fastq.0=_2.excerpt.fastq.$$i) \
			--out $(@:.sam.0=.sam.$$i) \
			-- -t 8 ; \
	done

# BWA-MEM unpaired overall
%.bwa.unp.sam: %.fastq
	python ../ts.py \
		$(BWA_ARGS) $(ARGS) \
		--U $< \
		--vanilla-output $@ \
		-- -t 8

# BWA-MEM unpaired samples
%.bwa.unp.sam.0: %.fastq.0
	for i in `seq 0 9`; do \
		python ../ts.py \
			$(BWA_ARGS) $(ARGS) $(ARGS_SAMP) \
			--U $(<:.0=.$$i) \
			--out $(@:.sam=.sam.$$i) \
			-- -t 8 ; \
	done

# #######
# Bowtie2
# #######

# Bowtie2 paired-end overall
%.bt2.pair.sam: %.fastq
	python ../ts.py \
		$(BT2_ARGS) $(ARGS) \
		--m1 $< \
		--m2 $(<:_1.excerpt.fastq=_2.excerpt.fastq) \
		--vanilla-output $@ \
		-- -p 8

# Bowtie2 paired-end samples
%.bt2.pair.sam.0: %.fastq.0
	for i in `seq 0 9`; do \
		python ../ts.py \
			$(BT2_ARGS) $(ARGS) $(ARGS_SAMP) \
			--m1 $(<:.0=.$$i) \
			--m2 $(<:_1.excerpt.fastq.0=_2.excerpt.fastq.$$i) \
			--out $(@:.sam=) \
			-- -p 8 ; \
	done

# Bowtie2 unpaired overall
%.bt2.unp.sam: %.fastq
	python ../ts.py \
		$(BT2_ARGS) $(ARGS) \
		--U $< \
		--vanilla-output $@ \
		-- -p 8

# Bowtie2 unpaired samples
%.bt2.unp.sam.0: %.fastq.0
	for i in `seq 0 9`; do \
		python ../ts.py \
			$(BT2_ARGS) $(ARGS) $(ARGS_SAMP) \
			--U $(<:.0=.$$i) \
			--out $(@:.sam=.sam.$$i) \
			-- -p 8 ; \
	done

# #######
# SNAP
# #######

# SNAP paired-end overall
%.snap.pair.sam: %.fastq
	python ../ts.py \
		$(SNAP_ARGS) $(ARGS) \
		--m1 $< \
		--m2 $(<:_1.excerpt.fastq=_2.excerpt.fastq) \
		--vanilla-output $@ \
		-- -t 8

# SNAP paired-end samples
%.snap.pair.sam.0: %.fastq.0
	for i in `seq 0 9`; do \
		python ../ts.py \
			$(SNAP_ARGS) $(ARGS) $(ARGS_SAMP) \
			--m1 $(<:.0=.$$i) \
			--m2 $(<:_1.excerpt.fastq.0=_2.excerpt.fastq.$$i) \
			--out $(@:.sam=.sam.$$i) \
			-- -t 8 ; \
	done

# SNAP unpaired overall
%.snap.unp.sam: %.fastq
	python ../ts.py \
		$(SNAP_ARGS) $(ARGS) \
		--U $< \
		--vanilla-output $@ \
		-- -t 8

# SNAP unpaired samples
%.snap.unp.sam.0: %.fastq.0
	for i in `seq 0 9`; do \
		python ../ts.py \
			$(SNAP_ARGS) $(ARGS) $(ARGS_SAMP) \
			--U $(<:.0=.$$i) \
			--out $(@:.sam=.sam.$$i) \
			-- -t 8 ; \
	done

# ###############
# Extra sensitive
# ###############

%.bt2vs.unp.sam: %.fastq
	python ../ts.py \
		$(BT2_ARGS) $(ARGS) \
		--U $< \
		--vanilla-output $@ \
		-- -p 8 --very-sensitive

%.bt2vsl.unp.sam: %.fastq
	python ../ts.py \
		$(BT2_ARGS) $(ARGS) \
		--U $< \
		--vanilla-output $@ \
		-- -p 8 --very-sensitive-local

%.fastq.0: %.fastq
	python sample_fastq.py -n 100000 -s 10 $< $(@:.0=) -r 43547
	python sample_fastq.py -n 100000 -s 10 $(<:_1.excerpt.fastq=_2.excerpt.fastq) $(@:_1.excerpt.fastq.0=_2.excerpt.fastq) -r 43547

.PHONY: clean
clean:
	rm -f snap_*.sam bwamem_*.sam bt2_*.sam
	rm -rf bt2_unpaired bwamem_paired bwamem_unpaired
