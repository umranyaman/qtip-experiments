# qtip-experiments

Scripts for driving all experiments described in the `qtip` manuscript.

### Clone repos

Clone the [`qtip`](https://github.com/BenLangmead/qtip) and [`qtip-experiments`](https://github.com/BenLangmead/qtip-experiments) repos.  `qtip` contains the Qtip software and scripts for building `qtip`-compatible versions of the relevant aligners.  `qtip-experiments` contains everything else needed to run the experiments described in the manuscript.

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

Make the software required for read alignment.  The read aligners are patched to be `qtip`-compatible:

```
make -C qtip/software/bowtie2
make -C qtip/software/bwa
make -C qtip/software/snap
```

Make the software required for read simulation:

```
make -f Makefile.src_linux -C qtip-experiments/software/art
make -C qtip-experiments/software/mason
make -C qtip-experiments/software/wgsim
```

Make the software required for evaluating variant calls:

```
make -C qtip-experiments/software/vcflib
```

MARCC note: making mason doesn’t work with `vtune` module loaded.

### Obtain reference genomes and build indexes

Reference genomes (`hg19`, `hg38`, `mm10` and `zm_AGPv4`) and indexes (`bowtie2`, `bwa mem` and `snap`) are used in these experiments.

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

### Obtain CHM1 assemblies and build combined references

```
pushd qtip-experiments/experiments/refs
# this may be slow, as it needs to download two human assemblies
sh get_assemblytics.sh
popd
```

### Alignment error experiments

For Supplementary Note 2.

```
pushd qtip-experiments/experiments/alignment_err
# obtain contaminant genomes
sh get_small_refs.sh
# simulate reads, compose target/foreign mixtures, align, summarize
make corstats
# tabulate
python tabulate.py
# analyze table
python evaluate.py > results.csv
```

### `simulated_reads` experiments

#### Simulate reads

```
pushd qtip-experiments/experiments/simulated_reads
python marcc_reads.py wet
# substitute "dry" for "wet" to just print the job-submission commands
popd
```

Many jobs are submitted here.  The longest job takes about 2 hours.

The script was written for the MARCC cluster at JHU; you might have to tweak for your cluster.

#### Run `qtip`

Scripts are described in more detail in the `qtip-experiments/experiments/simulated_reads/README.md` file.

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

### Training data formula series

To gather data on many ways of setting the `--sim-function` and `--sim-factor` parameters, which determine how many tandem reads to simulate as a function of the number of input reads, run:

```
pushd qtip-experiments/experiments/simulated_reads
sh train_series.sh --wet
```

Many jobs are submitted here.

The script was written for the MARCC cluster at JHU; you might have to tweak for your cluster.

### `platinum` experiments

#### Download reads

```
pushd qtip-experiments/experiments/platinum
sh get.sh
```

#### Download platinum variants and high confidence regions

```
sh get_gold_vcf.sh
```

#### Analyze reads

```
sh fraglen.sh wet
# previous is optional; useful for learning frag length dist
# pass "dry" instead to just print batch commands
sh align_full.sh wet
# pass "dry" instead to just print batch commands
# very time- and resource-intensive
```

#### Call variants

Calls variants:
* for all MAPQ thresholds
* for both original and qtip-predicted MAPQs
* for chromosomes 1-22 and X

```
sh sambamba_sort.sh
sh sbatch_fb.sh wet
# pass "dry" instead to just print batch commands
# very time- and resource-intensive
```

#### Summarize

Makes ROCs, and calculates F-scores based on best MAPQ and QUAL thresholds for both original and Qtip-predicted MAPQs.

```
sh sbatch_vcfroc.sh wet
python table_f.py > ERR194147.csv
popd
```
