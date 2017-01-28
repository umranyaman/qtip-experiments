Real data
---------

The `qtip` manuscript demonstrates how `qtip`-predicted mapping qualities improve performance when analyzing real data.

### Human

We use real data to demonstrate `qtip` performance.  Primarily we use two
human samples from the [1000 Genomes project](http://www.1000genomes.org):

* [ERR050082](http://www.ebi.ac.uk/ena/data/view/ERR050082)
    * 42,245,074 paired-end 100 x 100 nt reads
    * NA06985, CEPH, female
* [ERR050083](http://www.ebi.ac.uk/ena/data/view/ERR050083)
    * 66,391,067 paired-end 100 x 100 nt reads
    * NA11881, CEPH, male

Would be nice to know the sex here.

### Maize

(Future work)

Might also experiment with real data from a more repetitive genome, like maize.  If so, the links below might be relevant.

* http://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1592709
* ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/009/SRR1592709/SRR1592709_1.fastq.gz
* ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR159/009/SRR1592709/SRR1592709_2.fastq.gz

### Steps

Steps:

* Download FASTQ reads (`get_real_reads.sh`)
* Align FASTQ reads with various tools, producing SAM files (`sbatch_align.sh` / `Makefile`)
* Generate table of `qtip` overhead measurements (`perf_tabulate.py > perf.csv`)
* Create new reads annotated with "correctish" information using `sbatch_correctish.sh`
* Align the new, correctish-annotated reads using `sbatch_align_correctish.sh`
    * Once completed, the results are in the `new_*.sam` subdirectories
    * Eventually need to gather these up, like we do with `gather.py`

From directory containing the `qtip` and `qtip-experiments` repo clones:

```
pushd qtip-experiments/experiments/real_data
sh get_real_reads.sh  # might want to submit to your DRM
sh sbatch_align.sh
# copy and paste all the alignment jobs to submit them
# ...when those are done, proceed
sh sbatch_multialign.sh
# copy and paste all the alignment jobs to submit them
python perf_tabulate.py > perf.csv
python overall_tabulate.py
popd
```
