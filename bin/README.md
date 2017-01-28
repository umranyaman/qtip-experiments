These scripts are used from various other places in the repo, including from the `Makefile`s under `data/simulated_data` and from the SNAP test scripts under `software/snap`.

* `art_convert.py`: Convert Art-formatted FASTQ files to the augmented `wgsim`-like formatting used by the simulation scripts
* `mason_convert.py`: Same as above but for Mason-formatted FASTQ files
* `fastq_interleave.py`: Interleave two paired-end FASTQ files.  Sometimes useful for tools like BWA-MEM and SNAP that take interleaved FASTQ.
