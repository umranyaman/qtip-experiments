# qsim-experiments

The manuscript describing `qsim` describes several experiments.  This repo contains scripts for driving all these experiments.

## Simulated reads

Several of the MAPQ estimation accuracy results come from simulation experiments.  These experiments span a range of read lengths, aligners, simulators, reference genomes, etc.  All these experiments are driven by the `gmake` `Makefile`s in the `experiments/simulated_reads` subdirectory here.

To minimize overall time required, we ran many of these jobs in parallel on [MARCC](https://www.marcc.jhu.edu).  The scripts that drive the parallel jobs are described in the `experiments/simulated_reads/README.md` file.

## Real reads

We use real data to explore some aspects of qsim performance.  See the `experiments/real_data/README.md` file for details.

## Reference genomes and indexes

To run these experiments, the necessary reference genomes (`hg19`, `mm10` and `zm_AGPv3`) and indexes (`bowtie2`, `bwa mem` and `snap`) must be available.

The `experiments/refs/get_refs.sh` script downloads all the relevant reference genomes.  It also runs the `/experiments/refs/remove_short.py` script on some of them, which removes some short contigs that cause issues for certain read simulators.

The `experiments/refs/submit_index.sh` script submits nine index-building jobs, one for each aligner/genome combination.  The script was written for the HHPC cluster at JHU.

## Environment variables

For me, the path to my clone of the `qsim` GitHub repo is at `/scratch/users/blangme2@jhu.edu/git/mapq`, my genome indexes are in: `/scratch/groups/blangme2/indexes`, and my reference genomes are in: `/scratch/groups/blangme2/references`.  So I set up my `TS_*` environment variables like this:

```
TS_HOME=/scratch/users/blangme2@jhu.edu/git/mapq
TS_INDEXES=/scratch/groups/blangme2/indexes
TS_REFS=/scratch/groups/blangme2/references
```
