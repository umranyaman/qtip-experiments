TS_BWA_HOME=$(QSIM_EXPERIMENTS_HOME)/software/bwa
TS_SNAP_HOME=$(QSIM_EXPERIMENTS_HOME)/software/snap
TS_BT2_HOME=$(QSIM_EXPERIMENTS_HOME)/software/bowtie2
INDEXES=$(QSIM_EXPERIMENTS_HOME)/experiments/refs
PYTHON=python -O

BT2_QSIM_ARGS=--aligner bowtie2 --bt2-exe $(TS_BT2_HOME)/bowtie2 --index $(INDEXES)/hg19.fa
BT2_ARGS=
BT2_PAIR_ARGS=-X 750
BWA_QSIM_ARGS=--aligner bwa-mem --bwa-exe $(TS_BWA_HOME)/bwa --index $(INDEXES)/hg19.fa
BWA_ARGS=
BWA_PAIR_ARGS=
SNAP_QSIM_ARGS=--aligner snap --snap-exe $(TS_SNAP_HOME)/snap-aligner --index $(INDEXES)/hg19.fa.snap
SNAP_ARGS=
SNAP_PAIR_ARGS=

NTHREADS=1
NTHREADS_VS=16

$(shell mkdir -p temp)
ARGS=--ref $(INDEXES)/hg19.fa --write-orig-mapq --write-precise-mapq --seed 346 --temp-directory temp
EXTRA_ARGS=

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
	    --name bt2vs bt2vsl bwamem-qsim bt2-qsim snap-qsim \
	    --tier 1 1 2 2 2 > $@

%.pair.csv: %.bt2vs.pair.sam %.bt2vsl.pair.sam %.bwa.pair.sam %.bt2.pair.sam %.snap.pair.sam
	pypy multi_aligner.py \
	    --sam $^ \
	    --prefix $(<:bt2vs.pair.sam=) \
	    --suffix ".pair" \
	    --name bt2vs bt2vsl bwamem-qsim bt2-qsim snap-qsim \
	    --tier 1 1 2 2 2 > $@

# #######
# BWA-MEM
# #######

%.bwa.unp.sam: %.fastq
	$(PYTHON) $(QSIM_HOME)/src/qsim \
		$(BWA_QSIM_ARGS) $(ARGS) $(EXTRA_ARGS) \
		--U $< \
		--vanilla-output $@ \
		-- -t $(NTHREADS) $(BWA_ARGS)

%.bwa.pair.sam: %.fastq
	$(PYTHON) $(QSIM_HOME)/src/qsim \
		$(BWA_QSIM_ARGS) $(ARGS) $(EXTRA_ARGS) \
		--m1 $< \
		--m2 $(<:_1.fastq=_2.fastq) \
		--vanilla-output $@ \
		-- -t $(NTHREADS) $(BWA_ARGS) $(BWA_PAIR_ARGS)

# #######
# Bowtie2
# #######

%.bt2.unp.sam: %.fastq
	$(PYTHON) $(QSIM_HOME)/src/qsim \
		$(BT2_QSIM_ARGS) $(ARGS) $(EXTRA_ARGS) \
		--U $< \
		--vanilla-output $@ \
		-- -p $(NTHREADS) $(BT2_ARGS)

%.bt2.pair.sam: %.fastq
	$(PYTHON) $(QSIM_HOME)/src/qsim \
		$(BT2_QSIM_ARGS) $(ARGS) $(EXTRA_ARGS) \
		--m1 $< \
		--m2 $(<:_1.fastq=_2.fastq) \
		--vanilla-output $@ \
		-- -p $(NTHREADS) $(BT2_ARGS) $(BT2_PAIR_ARGS)

# #######
# SNAP
# #######

%.snap.unp.sam: %.fastq
	$(PYTHON) $(QSIM_HOME)/src/qsim \
		$(SNAP_QSIM_ARGS) $(ARGS) $(EXTRA_ARGS) \
		--U $< \
		--vanilla-output $@ \
		-- -t $(NTHREADS) $(SNAP_ARGS)

%.snap.pair.sam: %.fastq
	$(PYTHON) $(QSIM_HOME)/src/qsim \
		$(SNAP_QSIM_ARGS) $(ARGS) $(EXTRA_ARGS) \
		--m1 $< \
		--m2 $(<:_1.fastq=_2.fastq) \
		--vanilla-output $@ \
		-- -t $(NTHREADS) $(SNAP_ARGS) $(SNAP_PAIR_ARGS)

# ###############
# Extra sensitive
# ###############

%.bt2vs.unp.sam: %.fastq
	$(TS_BT2_HOME)/bowtie2 \
		-x $(TS_INDEXES)/hg19.fa \
		-p $(NTHREADS_VS) \
		-U $< \
		--very-sensitive \
		$(BT2_ARGS) \
		-S $@

%.bt2vs.pair.sam: %.fastq
	$(TS_BT2_HOME)/bowtie2 \
		-x $(TS_INDEXES)/hg19.fa \
		-p $(NTHREADS_VS) \
		-1 $< -2 $(<:_1.fastq=_2.fastq) \
		--very-sensitive \
		$(BT2_ARGS) \
		$(BT2_PAIR_ARGS) \
		-S $@

# #####################
# Extra sensitive local
# #####################

%.bt2vsl.unp.sam: %.fastq
	$(TS_BT2_HOME)/bowtie2 \
		-x $(TS_INDEXES)/hg19.fa \
		-p $(NTHREADS_VS) \
		-U $< \
		--very-sensitive-local \
		$(BT2_ARGS) \
		-S $@

%.bt2vsl.pair.sam: %.fastq
	$(TS_BT2_HOME)/bowtie2 \
		-x $(TS_INDEXES)/hg19.fa \
		-p $(NTHREADS_VS) \
		-1 $< -2 $(<:_1.fastq=_2.fastq) \
		--very-sensitive-local \
		$(BT2_ARGS) \
		$(BT2_PAIR_ARGS) \
		-S $@
