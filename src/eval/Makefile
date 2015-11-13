TS_BWA_HOME=$(TS_HOME)/software/bwa
TS_SNAP_HOME=$(TS_HOME)/software/snap
TS_BT2_HOME=$(TS_HOME)/software/bowtie2

BT2_ARGS=--aligner bowtie2 --bt2-exe $(TS_BT2_HOME)/bowtie2-git/bowtie2 --index $(TS_INDEXES)/hg19.fa
BWA_ARGS=--aligner bwa-mem --bwa-exe $(TS_BWA_HOME)/bwa --index $(TS_INDEXES)/hg19.fa
SNAP_ARGS=--aligner snap --snap-exe $(TS_SNAP_HOME)/snap/snap-aligner --index $(TS_INDEXES)/hg19.fa.snap

ARGS=--ref $(TS_REFS)/hg19.fa --write-orig-mapq --write-precise-mapq --seed 346
ARGS_SAMP=--keep-intermediates --predict-for-training

NTHREADS=16

.PHONY: all
all: ERR050082_1.unp.csv ERR050082_1.pair.csv \
     ERR050083_1.unp.csv ERR050083_1.pair.csv

%.unp.csv: %.bt2vs.unp.sam %.bt2vsl.unp.sam %.bwa.unp.sam %.bt2.unp.sam %.snap.unp.sam
	python multi_aligner.py \
	    --sam $^ \
	    --prefix $(<:.bwa.unp.sam=) \
	    --suffix ".unp" \
	    --name bt2vs bt2vsl bwamem bt2 snap \
	    --tier 1 1 2 2 2 > $@

%.pair.csv: %.bt2vs.pair.sam %.bt2vsl.pair.sam %.bwa.pair.sam %.bt2.pair.sam %.snap.pair.sam
	python multi_aligner.py \
	    --sam $^ \
	    --prefix $(<:.bwa.pair.sam=) \
	    --suffix ".pair" \
	    --name bt2vs bt2vsl bwamem bt2 snap \
	    --tier 1 1 2 2 2 > $@

# #######
# BWA-MEM
# #######

# BWA-MEM paired-end overall
%.bwa.pair.sam: %.fastq
	python ../ts.py \
		$(BWA_ARGS) $(ARGS) \
		--m1 $< \
		--m2 $(<:_1.fastq=_2.fastq) \
		--vanilla-output $@ \
		-- -t $(NTHREADS)

# BWA-MEM paired-end samples
%.bwa.pair.sam.0: %.fastq.0
	for i in `seq 0 9`; do \
		python ../ts.py \
			$(BWA_ARGS) $(ARGS) $(ARGS_SAMP) \
			--m1 $(<:.0=.$$i) \
			--m2 $(<:_1.fastq.0=_2.fastq.$$i) \
			--out $(@:.sam.0=.sam.$$i) \
			-- -t $(NTHREADS) ; \
	done

# BWA-MEM unpaired overall
%.bwa.unp.sam: %.fastq
	python ../ts.py \
		$(BWA_ARGS) $(ARGS) \
		--U $< \
		--vanilla-output $@ \
		-- -t $(NTHREADS)

# BWA-MEM unpaired samples
%.bwa.unp.sam.0: %.fastq.0
	for i in `seq 0 9`; do \
		python ../ts.py \
			$(BWA_ARGS) $(ARGS) $(ARGS_SAMP) \
			--U $(<:.0=.$$i) \
			--out $(@:.sam=.sam.$$i) \
			-- -t $(NTHREADS) ; \
	done

# #######
# Bowtie2
# #######

# Bowtie2 paired-end overall
%.bt2.pair.sam: %.fastq
	python ../ts.py \
		$(BT2_ARGS) $(ARGS) \
		--m1 $< \
		--m2 $(<:_1.fastq=_2.fastq) \
		--vanilla-output $@ \
		-- -p $(NTHREADS)

# Bowtie2 paired-end samples
%.bt2.pair.sam.0: %.fastq.0
	for i in `seq 0 9`; do \
		python ../ts.py \
			$(BT2_ARGS) $(ARGS) $(ARGS_SAMP) \
			--m1 $(<:.0=.$$i) \
			--m2 $(<:_1.fastq.0=_2.fastq.$$i) \
			--out $(@:.sam=) \
			-- -p $(NTHREADS) ; \
	done

# Bowtie2 unpaired overall
%.bt2.unp.sam: %.fastq
	python ../ts.py \
		$(BT2_ARGS) $(ARGS) \
		--U $< \
		--vanilla-output $@ \
		-- -p $(NTHREADS)

# Bowtie2 unpaired samples
%.bt2.unp.sam.0: %.fastq.0
	for i in `seq 0 9`; do \
		python ../ts.py \
			$(BT2_ARGS) $(ARGS) $(ARGS_SAMP) \
			--U $(<:.0=.$$i) \
			--out $(@:.sam=.sam.$$i) \
			-- -p $(NTHREADS) ; \
	done

# #######
# SNAP
# #######

# SNAP paired-end overall
%.snap.pair.sam: %.fastq
	python ../ts.py \
		$(SNAP_ARGS) $(ARGS) \
		--m1 $< \
		--m2 $(<:_1.fastq=_2.fastq) \
		--vanilla-output $@ \
		-- -t $(NTHREADS)

# SNAP paired-end samples
%.snap.pair.sam.0: %.fastq.0
	for i in `seq 0 9`; do \
		python ../ts.py \
			$(SNAP_ARGS) $(ARGS) $(ARGS_SAMP) \
			--m1 $(<:.0=.$$i) \
			--m2 $(<:_1.fastq.0=_2.fastq.$$i) \
			--out $(@:.sam=.sam.$$i) \
			-- -t $(NTHREADS) ; \
	done

# SNAP unpaired overall
%.snap.unp.sam: %.fastq
	python ../ts.py \
		$(SNAP_ARGS) $(ARGS) \
		--U $< \
		--vanilla-output $@ \
		-- -t $(NTHREADS)

# SNAP unpaired samples
%.snap.unp.sam.0: %.fastq.0
	for i in `seq 0 9`; do \
		python ../ts.py \
			$(SNAP_ARGS) $(ARGS) $(ARGS_SAMP) \
			--U $(<:.0=.$$i) \
			--out $(@:.sam=.sam.$$i) \
			-- -t $(NTHREADS) ; \
	done

# ###############
# Extra sensitive
# ###############

%.bt2vs.unp.sam: %.fastq
	python ../ts.py \
		$(BT2_ARGS) $(ARGS) \
		--U $< \
		--vanilla-output $@ \
		-- -p $(NTHREADS) --very-sensitive

%.bt2vs.pair.sam: %.fastq
	python ../ts.py \
		$(BT2_ARGS) $(ARGS) \
		--m1 $< \
		--m2 $(<:_1.fastq=_2.fastq) \
		--vanilla-output $@ \
		-- -p $(NTHREADS) --very-sensitive

%.bt2vsl.unp.sam: %.fastq
	python ../ts.py \
		$(BT2_ARGS) $(ARGS) \
		--U $< \
		--vanilla-output $@ \
		-- -p $(NTHREADS) --very-sensitive-local

%.bt2vsl.pair.sam: %.fastq
	python ../ts.py \
		$(BT2_ARGS) $(ARGS) \
		--m1 $< \
		--m2 $(<:_1.fastq=_2.fastq) \
		--vanilla-output $@ \
		-- -p $(NTHREADS) --very-sensitive-local

%.fastq.0: %.fastq
	python sample_fastq.py -n 100000 -s 10 $< $(@:.0=) -r 43547
	python sample_fastq.py -n 100000 -s 10 $(<:_1.fastq=_2.fastq) $(@:_1.fastq.0=_2.fastq) -r 43547

.PHONY: clean
clean:
	rm -f *.unp.sam *.pair.sam
	rm -f *_roc.csv *_summ.csv
	rm -f *.unp.csv *.pair.csv
	rm -rf bt2_unpaired bwamem_paired bwamem_unpaired
