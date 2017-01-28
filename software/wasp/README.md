WASP
====

[WASP] uses [Snakemake] and a variety of Python libraries.
To run WASP, you will first need to switch into a python3 environment.
On my machine, I do this with `source activate py3k` for example.

Snakemake uses a YAML configuration file.
When you first clone WASP, the configuration file has someone else's directories in it.
So the `run_snakemake.sh` script won't work until you've corrected that.

I tried to get Snakemake to work with the example data that comes with WASP.

[WASP]: https://github.com/bmvdgeijn/WASP
[Snakemake]: https://bitbucket.org/snakemake/

Mapping
-------

At the end of the day, the part of WASP we're most concerned with is the component that removes bias by re-mapping reads to the genome with variants inserted.
This software is in the `mapping` subdirectory of the repo.
Here are some of the relevant scripts and some notes about them

#### `find_intersecting_snps.py`

Depends on `tables` and `snptables` modules?  Are these PyTables?
`SNPTable` class is in `snptable.py`.

Needs input to be sorted.
It loads in SNPs in chromosome-by-chromosome batches.

The crux of it is in `process_single_read` and `process_paired_read`.
Calls `generate_reads` or `generate_haplo_reads` as needed.

Basically it:
* Sorts input SAM/BAM into sorted BAM, if needed

If I want to use it with some human RNA sequencing data (or any other kind of data) and the human GRCh38 genome, then I need to do the following things:

* Get variants for some population, with GRCh38 coordinates
* Arrange them as needed in `--snp_dir`
    * Named like `chr<#>.snps.txt.gz`
    * 3 columns: position, RefAllele, AltAllele
* Give it a BAM file with `--bam_filename`
* Set `--is_paired_end` if appropriate

#### `filter_remapped_reads.py`

Runs after `find_intersecting_snps.py`.