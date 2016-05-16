# qsim-experiments

Scripts for driving all experiments described in the `qsim` manuscript.

### Clone repos

Clone the [`qsim`](https://github.com/BenLangmead/qsim) and [`qsim-experiments`](https://github.com/BenLangmead/qsim-experiments) repos.  The `qsim` repo contains the `qsim` software and scripts for building `qsim`-compatible versions of the relevant read aligners.  `qsim-experiments` contains everything else needed to run the experiments described in the manuscript.

```
cd $US/git
git clone git@github.com:BenLangmead/qsim.git
git clone git@github.com:BenLangmead/qsim-experiments.git
```

Substitute an appropriate directory for `$US/git` above.  `$US/git` is the user scratch directory on our MARCC cluster.

### Set up environment

Assuming you're still in the directory from which you did the `git clone`s:

```
export QSIM_HOME=`pwd`/qsim
export QSIM_EXPERIMENTS_HOME=`pwd`/qsim-experiments
```

Consider adding corresponding commands to `.bashrc` or the like.

### Build `qsim`

```
make -C $QSIM_HOME/src
$QSIM_HOME/src/qsim --version
```

`make` is required because, although most of `qsim` is in Python, some helper programs are C++.  You'll need a C++ compiler for this.

### Build software in `qsim-experiments`

Make the required software for read alignment and read simulation.  The read aligners are patched to be `qsim`-compatible.

```
make -f Makefile.src_linux -C qsim-experiments/software/art
make -C qsim-experiments/software/mason
make -C qsim-experiments/software/wgsim
make -C qsim-experiments/software/bowtie2
make -C qsim-experiments/software/bwa
make -C qsim-experiments/software/snap
```

MARCC note: making mason doesn’t work with `vtune` module loaded.

### Obtain reference genomes and build indexes

Reference genomes (`hg19`, `mm10` and `zm_AGPv3`) and indexes (`bowtie2`, `bwa mem` and `snap`) are used in these experiments.

```
pushd qsim-experiments/experiments/refs
sh get_refs.sh
```

`get_refs.sh` both downloads the relevant reference genomes and runs the `qsim-experiments/experiments/refs/remove_short.py` script on some of them, removing short contigs that cause issues for certain read simulators.

The `submit_index.sh` script submits nine index-building jobs, one for each aligner/genome combination.  The script was written for the MARCC cluster at JHU; you might have to tweak for your cluster.

```
sh index_marcc.sh
# copy and paste the printed sbatch commands to actually submit them
# these jobs will take a few hours
popd
```

### Simulate reads in `qsim-experiments`

```
pushd qsim-experiments/experiments/simulated_reads
python marcc_reads.py wet
# substitute "dry" for "wet" to just print the job-submission commands
popd
```

Many jobs are submitted here.  The longest job takes about 2 hours.

The script was written for the MARCC cluster at JHU; you might have to tweak for your cluster.

### Run `qsim` on simulated datasets in `qsim-experiments`

These scripts are described in more detail in the `qsim-experiments/experiments/simulated_reads/README.md` file.

```
pushd qsim-experiments/experiments/simulated_reads
python marcc_out.py wet
# substitute "dry" for "wet" to just print the job-submission commands
popd
```

Many jobs are submitted here.  All told, this takes about 4 hours for me on MARCC.

The script was written for the MARCC cluster at JHU; you might have to tweak for your cluster.

### Run qsim on real datasets in “qsim-experiments”

See `qsim-experiments/experiments/real_data/README.md` file for more details.

Note that the user has to issue the commands to the distributed resource manager.  The `sbatch_all.sh` script generates the scripts for the jobs and prints out appropriate submission commands assuming the SLURM DRM, as exists on MARCC.  But you may have to modify for your cluster.

```
pushd qsim-experiments/experiments/real_data
sh get_real_reads.sh  # might want to submit to your DRM
sh sbatch_all.sh
# copy and paste all the alignment jobs to submit them
# ...when those are done, copy and paste all the multi_aligner.py jobs
popd
```
