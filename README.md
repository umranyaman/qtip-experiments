# qsim-experiments

The manuscript describing `qsim` describes several experiments.  This repo contains scripts for driving all these experiments.

## Simulated reads

Several of the MAPQ estimation accuracy results come from simulation experiments.  These experiments span a range of read lengths, aligners, simulators, reference genomes, etc.  All these experiments are driven by the `gmake` `Makefile`s in the `experiments/simulated_reads` subdirectory here.

To minimize overall time required, we ran many of these jobs in parallel on [MARCC](https://www.marcc.jhu.edu).  The scripts use to drive these experiments are described in somewhat more detail in the `experiments/simulated_reads/README.md` file.
