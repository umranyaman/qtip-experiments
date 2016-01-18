TS_BWA_HOME=$(TS_HOME)/software/bwa
TS_SNAP_HOME=$(TS_HOME)/software/snap
TS_BT2_HOME=$(TS_HOME)/software/bowtie2
PYTHON=python -O

BT2_ARGS=--aligner bowtie2 --bt2-exe $(TS_BT2_HOME)/bowtie2-git/bowtie2 --index $(TS_INDEXES)/hg19.fa
BWA_ARGS=--aligner bwa-mem --bwa-exe $(TS_BWA_HOME)/bwa --index $(TS_INDEXES)/hg19.fa
SNAP_ARGS=--aligner snap --snap-exe $(TS_SNAP_HOME)/snap/snap-aligner --index $(TS_INDEXES)/hg19.fa.snap

ARGS=--ref $(TS_REFS)/hg19.fa --write-orig-mapq --write-precise-mapq --seed 346

NTHREADS=1
NTHREADS_VS=16

#
# These CSV files contain the ultimate output from multi_aligner.py when
# analyzing SAMs from various aligners (bt2 --very-sensitive,
# bt2 --very-sensitive-local, bwa-mem, bt2 --sensitive, and SNAP).
#

.PHONY: all
all: ERR050082_1.unp.csv ERR050082_1.pair.csv \
     ERR050083_1.unp.csv ERR050083_1.pair.csv

%.unp.csv: %.bt2vs.unp.sam %.bt2vsl.unp.sam %.bwa.unp.sam %.bt2.unp.sam %.snap.unp.sam
	pypy multi_aligner.py \
	    --sam $^ \
	    --prefix $(<:bt2vs.unp.sam=) \
	    --suffix ".unp" \
	    --name bt2vs bt2vsl bwamem bt2 snap \
	    --tier 1 1 2 2 2 > $@

%.pair.csv: %.bt2vs.pair.sam %.bt2vsl.pair.sam %.bwa.pair.sam %.bt2.pair.sam %.snap.pair.sam
	pypy multi_aligner.py \
	    --sam $^ \
	    --prefix $(<:bt2vs.pair.sam=) \
	    --suffix ".pair" \
	    --name bt2vs bt2vsl bwamem bt2 snap \
	    --tier 1 1 2 2 2 > $@

# #######
# BWA-MEM
# #######

%.bwa.unp.sam: %.fastq
	$(PYTHON) ../ts.py \
		$(BWA_ARGS) $(ARGS) \
		--U $< \
		--vanilla-output $@ \
		-- -t $(NTHREADS)

%.bwa.pair.sam: %.fastq
	$(PYTHON) ../ts.py \
		$(BWA_ARGS) $(ARGS) \
		--m1 $< \
		--m2 $(<:_1.fastq=_2.fastq) \
		--vanilla-output $@ \
		-- -t $(NTHREADS)

#%.bwa.pair.sam.plain: %.fastq
#	$(TS_BWA_HOME)/bwa \
#		mem \
#		-t $(NTHREADS) \
#		$(TS_INDEXES)/hg19.fa \
#		$< $(<:_1.fastq=_2.fastq) \
#		> $@

#%.bwa.unp.sam.plain: %.fastq
#	$(TS_BWA_HOME)/bwa \
#		mem \
#		-t $(NTHREADS) \
#		$(TS_INDEXES)/hg19.fa \
#		$< \
#		> $@

# #######
# Bowtie2
# #######

%.bt2.unp.sam: %.fastq
	$(PYTHON) ../ts.py \
		$(BT2_ARGS) $(ARGS) \
		--U $< \
		--vanilla-output $@ \
		-- -p $(NTHREADS)

%.bt2.pair.sam: %.fastq
	$(PYTHON) ../ts.py \
		$(BT2_ARGS) $(ARGS) \
		--m1 $< \
		--m2 $(<:_1.fastq=_2.fastq) \
		--vanilla-output $@ \
		-- -p $(NTHREADS)

#%.bt2.pair.sam.plain: %.fastq
#	$(TS_BT2_HOME)/bowtie2 \
#		-x $(TS_INDEXES)/hg19.fa \
#		-p $(NTHREADS) \
#		-1 $< -2 $(<:_1.fastq=_2.fastq) \
#		-S $@

#%.bt2.unp.sam.plain: %.fastq
#	$(TS_BT2_HOME)/bowtie2 \
#		-x $(TS_INDEXES)/hg19.fa \
#		-p $(NTHREADS) \
#		-U $< \
#		-S $@

# #######
# SNAP
# #######

%.snap.unp.sam: %.fastq
	$(PYTHON) ../ts.py \
		$(SNAP_ARGS) $(ARGS) \
		--U $< \
		--vanilla-output $@ \
		-- -t $(NTHREADS)

%.snap.pair.sam: %.fastq
	$(PYTHON) ../ts.py \
		$(SNAP_ARGS) $(ARGS) \
		--m1 $< \
		--m2 $(<:_1.fastq=_2.fastq) \
		--vanilla-output $@ \
		-- -t $(NTHREADS)

#%.snap.pair.sam.plain: %.fastq
#	$(TS_SNAP_HOME)/snap/snap-aligner \
#		paired \
#		$(TS_INDEXES)/hg19.fa.snap \
#		-fastq $< $(<:_1.fastq=_2.fastq) \
#		-o -sam $@ \
#		-t $(NTHREADS) \
#		-=

#%.snap.unp.sam.plain: %.fastq
#	$(TS_SNAP_HOME)/snap/snap-aligner \
#		single \
#		$(TS_INDEXES)/hg19.fa.snap \
#		-fastq $< \
#		-o -sam $@ \
#		-t $(NTHREADS) \
#		-=

# ###############
# Extra sensitive
# ###############

#%.bt2vs.unp.sam: %.fastq
#	$(PYTHON) ../ts.py \
#		$(BT2_ARGS) $(ARGS) \
#		--U $< \
#		--vanilla-output $@ \
#		-- -p $(NTHREADS) --very-sensitive

#%.bt2vs.pair.sam: %.fastq
#	$(PYTHON) ../ts.py \
#		$(BT2_ARGS) $(ARGS) \
#		--m1 $< \
#		--m2 $(<:_1.fastq=_2.fastq) \
#		--vanilla-output $@ \
#		-- -p $(NTHREADS) --very-sensitive

%.bt2vs.unp.sam: %.fastq
	$(TS_BT2_HOME)/bowtie2 \
		-x $(TS_INDEXES)/hg19.fa \
		-p $(NTHREADS_VS) \
		-U $< \
		--very-sensitive \
		-S $@

%.bt2vs.pair.sam: %.fastq
	$(TS_BT2_HOME)/bowtie2 \
		-x $(TS_INDEXES)/hg19.fa \
		-p $(NTHREADS_VS) \
		-1 $< -2 $(<:_1.fastq=_2.fastq) \
		--very-sensitive \
		-S $@

# #####################
# Extra sensitive local
# #####################

#%.bt2vsl.unp.sam: %.fastq
#	$(PYTHON) ../ts.py \
#		$(BT2_ARGS) $(ARGS) \
#		--U $< \
#		--vanilla-output $@ \
#		-- -p $(NTHREADS) --very-sensitive-local

#%.bt2vsl.pair.sam: %.fastq
#	$(PYTHON) ../ts.py \
#		$(BT2_ARGS) $(ARGS) \
#		--m1 $< \
#		--m2 $(<:_1.fastq=_2.fastq) \
#		--vanilla-output $@ \
#		-- -p $(NTHREADS) --very-sensitive-local

%.bt2vsl.unp.sam: %.fastq
	$(TS_BT2_HOME)/bowtie2 \
		-x $(TS_INDEXES)/hg19.fa \
		-p $(NTHREADS_VS) \
		-U $< \
		--very-sensitive-local \
		-S $@

%.bt2vsl.pair.sam: %.fastq
	$(TS_BT2_HOME)/bowtie2 \
		-x $(TS_INDEXES)/hg19.fa \
		-p $(NTHREADS_VS) \
		-1 $< -2 $(<:_1.fastq=_2.fastq) \
		--very-sensitive-local \
		-S $@
