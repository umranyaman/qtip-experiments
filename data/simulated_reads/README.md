Simulated read datasets
-----------------------

Subdirectories contain sets of reads needed for various experiments:

* `ill_various_length` Illumina-like reads simulated by Mason.  Reads
  are in many fastq files with all reads of the same length in a
  single file.  Lengths simulated were: 50, 100, 150, 250, 500.  Both
  unpaired and paried-end reads simulated for each length.  hg19 is
  the reference simulated from.
* `simple`: like `ill_various_length` but with just 100 and 250 nt
  reads.  Useful for experiments where the aligner or the aligner
  parameters are being varied.
* `454_mixed_length`: one fastq file with 454 reads simulated from
  Mason.  Read lengths are drawn uniformly from [50, 500].  Unpaired.
  hg19 is the reference simulated from.
* `various_genomes`: like `simple` but where the genome is being varied.
* `various_simulators`
