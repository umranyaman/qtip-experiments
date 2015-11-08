TS_BWA_HOME=$(TS_HOME)/software/bwa
TS_SNAP_HOME=$(TS_HOME)/software/snap
TS_BT2_HOME=$(TS_HOME)/software/bowtie2

BT2_ARGS=--aligner bowtie2 --bt2-exe $(TS_BT2_HOME)/bowtie2-git/bowtie2 --index $(TS_INDEXES)/hg19.fa
BWA_ARGS=--aligner bwa-mem --bwa-exe $(TS_BWA_HOME)/bwa --index $(TS_INDEXES)/hg19.fa
SNAP_ARGS=--aligner snap --snap-exe $(TS_SNAP_HOME)/snap/snap-aligner --index $(TS_INDEXES)/hg19.fa.snap

ARGS=--ref $(TS_REFS)/hg19.fa --write-orig-mapq --write-precise-mapq --keep-intermediates --seed 346 --predict-for-training

.PHONY: all
all: even_tiered_unpaired.csv even_tiered_paired.csv

even_tiered_unpaired.csv: bt2vs_unpaired.sam bt2vsl_unpaired.sam bwamem_unpaired.sam bt2_unpaired.sam snap_unpaired.sam
	python multi_aligner.py \
	    --sam $^ \
	    --name bt2vs bt2vsl bwamem bt2 snap \
	    --tier 1 1 2 2 2 > $@

even_tiered_paired.csv: bwamem_paired.sam bt2_paired.sam snap_paired.sam
	python multi_aligner.py --sam $^ --tier 1 1 1 > $@

%.bwa.pair.sam: %.fastq.1
	for i in `seq 0 9`; do \
		python ../ts.py \
			$(BWA_ARGS) $(ARGS) \
			--m1 $(<:.1=.$$i) \
			--m2 $(<:_1.excerpt.fastq.1=_2.excerpt.fastq.$$i) \
			--out $(@:.sam=.sam.$$i) \
			-- -t 8 ; \
	done

%.bwa.unp.sam: %.fastq.1
	for i in `seq 0 9`; do \
		python ../ts.py \
			$(BWA_ARGS) $(ARGS) \
			--U $(<:.1=.$$i) \
			--out $(@:.sam=.sam.$$i) \
			-- -t 8 ; \
	done

%.bt2.pair.sam: %.fastq.1
	for i in `seq 0 9`; do \
		python ../ts.py \
			$(BT2_ARGS) $(ARGS) \
			--m1 $(<:.1=.$$i) \
			--m2 $(<:_1.excerpt.fastq.1=_2.excerpt.fastq.$$i) \
			--out $(@:.sam=) \
			-- -p 8 ; \
	done

%.bt2.unp.sam: %.fastq.1
	for i in `seq 0 9`; do \
		python ../ts.py \
			$(BT2_ARGS) $(ARGS) \
			--U $(<:.1=.$$i) \
			--out $(@:.sam=.sam.$$i) \
			-- -p 8 ; \
	done

# Extra sensitive
bt2vs_unpaired.sam: ERR050082_1.sm_excerpt.fastq
	python ../ts.py \
		$(BT2_ARGS) $(ARGS) \
		--U $< \
		--out $(@:.sam=) \
		-- -p 8

%.snap.unp.sam: %.fastq.1
	for i in `seq 0 9`; do \
		python ../ts.py \
			$(SNAP_ARGS) $(ARGS) \
			--U $(<:.1=.$$i) \
			--out $(@:.sam=.sam.$$i) \
			-- -t 8 ; \
	done

%.snap.pair.sam: %.fastq.1
	for i in `seq 0 9`; do \
		python ../ts.py \
			$(SNAP_ARGS) $(ARGS) \
			--m1 $(<:.1=.$$i) \
			--m2 $(<:_1.excerpt.fastq.1=_2.excerpt.fastq.$$i) \
			--out $(@:.sam=) \
			-- -t 8 ; \
	done

%.fastq.1: %.fastq
	python sample_fastq.py -n 100000 -s 10 $< $(@:.1=) -r 43547
	python sample_fastq.py -n 100000 -s 10 $(<:_1.excerpt.fastq=_2.excerpt.fastq) $(@:_1.excerpt.fastq.1=_2.excerpt.fastq) -r 43547

.PHONY: clean
clean:
	rm -f snap_*.sam bwamem_*.sam bt2_*.sam
	rm -rf bt2_unpaired bwamem_paired bwamem_unpaired
