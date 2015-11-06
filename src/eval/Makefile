TS_BWA_HOME=$(TS_HOME)/software/bwa
TS_SNAP_HOME=$(TS_HOME)/software/snap
TS_BT2_HOME=$(TS_HOME)/software/bowtie2

.PHONY: all
all: bwamem_paired.sam bt2_paired.sam snap_paired.sam \
     bwamem_unpaired.sam bt2_unpaired.sam snap_unpaired.sam

bwamem_paired.sam: ERR050082_1.sm_excerpt.fastq ERR050082_2.sm_excerpt.fastq
	python ../ts.py \
		--ref $(TS_REFS)/hg19.fa \
		--m1 ERR050082_1.sm_excerpt.fastq \
		--m2 ERR050082_2.sm_excerpt.fastq \
		--index $(TS_INDEXES)/hg19.fa \
		--seed 34 \
		--aligner bwa-mem \
		--bwa-exe $(TS_BWA_HOME)/bwa \
		--vanilla-output $@ \
		--write-orig-mapq \
		--write-precise-mapq

bwamem_unpaired.sam: ERR050082_1.sm_excerpt.fastq
	python ../ts.py \
		--ref $(TS_REFS)/hg19.fa \
		--m1 $< \
		--index $(TS_INDEXES)/hg19.fa \
		--seed 34 \
		--aligner bwa-mem \
		--bwa-exe $(TS_BWA_HOME)/bwa \
		--vanilla-output $@ \
		--write-orig-mapq \
		--write-precise-mapq

bt2_paired.sam: ERR050082_1.sm_excerpt.fastq ERR050082_2.sm_excerpt.fastq
	python ../ts.py \
		--ref $(TS_REFS)/hg19.fa \
		--m1 ERR050082_1.sm_excerpt.fastq \
		--m2 ERR050082_2.sm_excerpt.fastq \
		--index $(TS_INDEXES)/hg19.fa \
		--seed 34 \
		--aligner bowtie2 \
		--bt2-exe $(TS_BT2_HOME)/bowtie2 \
		--vanilla-output $@ \
		--write-orig-mapq \
		--write-precise-mapq

bt2_unpaired.sam: ERR050082_1.sm_excerpt.fastq,
	python ../ts.py \
		--ref $(TS_REFS)/hg19.fa \
		--m1 $< \
		--index $(TS_INDEXES)/hg19.fa \
		--seed 34 \
		--aligner bowtie2 \
		--bt2-exe $(TS_BT2_HOME)/bowtie2 \
		--vanilla-output $@ \
		--write-orig-mapq \
		--write-precise-mapq

snap_paired.sam: ERR050082_1.sm_excerpt.fastq ERR050082_2.sm_excerpt.fastq
	python ../ts.py \
		--ref $(TS_REFS)/hg19.fa \
		--m1 ERR050082_1.sm_excerpt.fastq \
		--m2 ERR050082_2.sm_excerpt.fastq \
		--index $(TS_INDEXES)/hg19.fa \
		--seed 34 \
		--aligner snap \
		--snap-exe $(TS_SNAP_HOME)/snap/snap-aligner \
		--vanilla-output $@ \
		--write-orig-mapq \
		--write-precise-mapq

snap_unpaired.sam: ERR050082_1.sm_excerpt.fastq
	python ../ts.py \
		--ref $(TS_REFS)/hg19.fa \
		--m1 $< \
		--index $(TS_INDEXES)/hg19.fa \
		--seed 34 \
		--aligner snap \
		--snap-exe $(TS_SNAP_HOME)/snap/snap-aligner \
		--vanilla-output $@ \
		--write-orig-mapq \
		--write-precise-mapq

%.sm_excerpt.fastq: %.excerpt.fastq
	head -n 100000 $< > $@

.PHONY: clean
clean:
	rm -f snap_*.sam bwamem_*.sam bt2_*.sam
