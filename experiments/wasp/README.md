## WASP experiments

WASP is essentially a scheme for rejecting alignments that seem likely to be incorrect.  It is proposed in the context of reducing bias in ASE & eQTL analysis.  But it is relevant to assessing mapping qualities too.  We ask: do mapping-quality thresholds provide an equally good or better way to decide which alignments to reject?  This could be assessed:

* In simulation, where a rejection decision is either correct or incorrect
* By assessing something downstream, like evenness of allele coverage or quality of var calls

### Real data experiments

* We create links to several relevant alignments files with Qtip predictions in `../real_data`
    * Must be from several aligners run on same input reads
* Genome-sort all of these alignments, in preparation for `find_intersecting_snps.py`
* `sbatch_wasp.sh` runs WASP's `find_intersecting_snps.py` script on each sorted BAM, giving FASTQ
* `sbatch_postprocess.sh` & `postprocess.py` do a few things:
    * Run the aligner on the WASP-generated FASTQs
    * Parse the SAM for the WASP-generated FASTQ
* ultimate output from `postprocess.py` is a bunch of records giving, for each input alignment, :
    * Original MAPQ
    * Qtip-predicted MAPQs
    * Number remapped reads aligned correctly
    * Number remapped reads aligned incorrectly
    * Number of times the previous four things co-occurred (for compression)

### Simulation experiments

* Alignments are copied from `../simulated_data/various_genomes`
* `sbatch_wasp.sh` runs WASP's `find_intersecting_snps.py` script on each, which generates FASTQ
* `sbatch_postprocess.sh` & `postprocess.py` look at the SAM and output; report is slightly different since the reads are simulated so we have truth info:
    * Original MAPQ
    * Qtip-predicted MAPQs
    * Number remapped reads aligned correctly
    * Number remapped reads aligned incorrectly
    * _Whether the alignment is correct_
    * Number of times the previous four things co-occurred (for compression)
