# qsim-experiments

Scripts for driving all experiments described in the `qsim` manuscript.

### Clone repos

```
cd $US/git  # substitute appropriately; this is "user scratch" on MARCC cluster
git clone git@github.com:BenLangmead/qsim.git
git clone git@github.com:BenLangmead/qsim-experiments.git
```

### Set up environment

```
export QSIM_HOME=`pwd`/qsim
export QSIM_EXPERIMENTS_HOME=`pwd`/qsim-experiments
```

Consider adding corresponding commands to `.bashrc` or the like.

### Build software in `qsim`

```
make -C $QSIM_HOME/src
$QSIM_HOME/src/qsim --version
```

`make` is required because, although most of `qsim` is in Python, some helper programs are C++.

### Build software in `qsim-experiments`

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

Reference genomes (`hg19`, `mm10` and `zm_AGPv3`) and indexes (`bowtie2`, `bwa mem` and `snap`) must be available for these experiments.

```
pushd qsim-experiments/experiments/refs
sh get_refs.sh
```

Note: `get_refs.sh` both downloads the relevant reference genomes and runs the `/experiments/refs/remove_short.py` script on some of them, removing short contigs that cause issues for certain read simulators.

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

See `qsim-experiments/experiments/real_data/README.md` file for details.

```
pushd qsim-experiments/experiments/real_data
sh get_real_reads.sh
sh sbatch_all.sh
popd
```
