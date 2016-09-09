# qsim-experiments

Scripts for driving all experiments described in the `qtip` manuscript.

### Clone repos

Clone the [`qtip`](https://github.com/BenLangmead/qtip) and [`qtip-experiments`](https://github.com/BenLangmead/qtip-experiments) repos.  The `qtip` repo contains the `qtip` software and scripts for building `qtip`-compatible versions of the relevant read aligners.  `qtip-experiments` contains everything else needed to run the experiments described in the manuscript.

```
cd $US/git
git clone git@github.com:BenLangmead/qtip.git
git clone git@github.com:BenLangmead/qtip-experiments.git
```

Substitute an appropriate directory for `$US/git` above.  `$US/git` is the user scratch directory on our MARCC cluster.

### Set up environment

Assuming you're still in the directory from which you did the `git clone`s:

```
export QTIP_HOME=`pwd`/qtip
export QTIP_EXPERIMENTS_HOME=`pwd`/qtip-experiments
```

Consider adding corresponding commands to `.bashrc` or the like.

### Build `qtip`

```
make -C $QTIP_HOME/src
$QTIP_HOME/src/qtip --version
```

Although most of `qtip` is in Python, some helper programs are C++.  You'll need a C++ compiler for this.

### Build software in `qtip-experiments`

Make the required software for read alignment and read simulation.  The read aligners are patched to be `qtip`-compatible.

```
make -f Makefile.src_linux -C qtip-experiments/software/art
make -C qtip-experiments/software/mason
make -C qtip-experiments/software/wgsim
make -C qtip-experiments/software/bowtie2
# For now: make -C qtip-experiments/software/bowtie2 -f Makefile.from_github
make -C qtip-experiments/software/bwa
make -C qtip-experiments/software/snap
```

MARCC note: making mason doesn’t work with `vtune` module loaded.

### Obtain reference genomes and build indexes

Reference genomes (`hg19`, `mm10` and `zm_AGPv3`) and indexes (`bowtie2`, `bwa mem` and `snap`) are used in these experiments.

```
pushd qtip-experiments/experiments/refs
sh get_refs.sh
```

`get_refs.sh` both downloads the relevant reference genomes and runs the `qtip-experiments/experiments/refs/remove_short.py` script on some of them, removing short contigs that cause issues for certain read simulators.

The `submit_index.sh` script submits nine index-building jobs, one for each aligner/genome combination.  The script was written for the MARCC cluster at JHU; you might have to tweak for your cluster.

```
sh index_marcc.sh
# copy and paste the printed sbatch commands to actually submit them
# these jobs will take a few hours
popd
```

### Simulate reads in `qtip-experiments`

```
pushd qtip-experiments/experiments/simulated_reads
python marcc_reads.py wet
# substitute "dry" for "wet" to just print the job-submission commands
popd
```

Many jobs are submitted here.  The longest job takes about 2 hours.

The script was written for the MARCC cluster at JHU; you might have to tweak for your cluster.

### Run `qtip` on simulated datasets in `qtip-experiments`

These scripts are described in more detail in the `qtip-experiments/experiments/simulated_reads/README.md` file.

```
pushd qtip-experiments/experiments/simulated_reads
python marcc_out.py wet
# substitute "dry" for "wet" to just print the job-submission commands
python gather.py --slurm
sbatch .gather.sh  # or edit as appropriate for your cluster
popd
```

Many jobs are submitted here.  All told, this takes about 4 hours for me on MARCC.

The script was written for the MARCC cluster at JHU; you might have to tweak for your cluster.

TODO: those steps do everything but the gathering of CID and CSED curves into a file that gets loaded into R.  Need at add that.

### Run `qtip` on real datasets in `qtip-experiments`

See `qtip-experiments/experiments/real_data/README.md` file for more details.

The user has to issue the commands to the distributed resource manager.  The `sbatch_align.sh` and `sbatch_multialign.sh` scripts generate the jobs and print appropriate submission commands without running them.

The `sbatch_*` scripts are intended for the SLURM DRM, which we use on the MARCC cluster.  You may have to modify the scripts for your cluster.

```
pushd qtip-experiments/experiments/real_data
sh get_real_reads.sh  # might want to submit to your DRM
sh sbatch_align.sh
# copy and paste all the alignment jobs to submit them
# ...when those are done, proceed
sh sbatch_multialign.sh
# copy and paste all the alignment jobs to submit them
popd
```
