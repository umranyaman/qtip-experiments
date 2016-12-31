## WASP experiments

WASP is essentially a scheme for rejecting alignments that seem likely to be incorrect.  It is proposed in the context of reducing bias in ASE & eQTL analysis.  But it is relevant to assessing mapping qualities too.  We ask: do mapping-quality thresholds provide an equally good or better way to decide which alignments to reject?  This could be assessed:

* In simulation, where a rejection decision is either correct or incorrect
* By assessing something downstream, like evenness of allele coverage or quality of var calls

### Real data experiments

* Alignments are copied from `../real_data`
* `sbatch_wasp.sh` runs WASP on each
* `sbatch_postprocess.sh` & `postprocess.py` look at the SAM and output; reports:
    * Original MAPQ
    * Qtip-predicted MAPQs
    * Number remapped reads aligned correctly
    * Number remapped reads aligned incorrectly
    * Number of times the previous four things co-occurred (for compression)

### Simulation experiments

* Alignments are copied from `../simulated_data/various_genomes`
* `sbatch_wasp.sh` runs WASP on each
* `sbatch_postprocess.sh` & `postprocess.py` look at the SAM and output; report is slightly different since the reads are simulated so we have truth info:
    * Original MAPQ
    * Qtip-predicted MAPQs
    * Number remapped reads aligned correctly
    * Number remapped reads aligned incorrectly
    * _Whether the alignment is correct_
    * Number of times the previous four things co-occurred (for compression)
