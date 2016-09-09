Contents
--------

* `marcc_reads.py` generates and optionally submits MARCC SLURM jobs for simulating reads
    * Parses the `Makefiles` in the subdirectories to figure out what jobs to submit
* `marcc_out.py` generates and optionally submits MARCC SLURM jobs for running Qtip+aligner
    * Parses the `Makefiles` in the subdirectories to figure out what jobs to submit
* `gather.py` parses contents of `.out` subdirectories

Helpful commands
----------------

To show all results in a particular subdirectory:

* `for i in `find . -name 'r*.out'` ; do echo $i ; cat $i/sample0.03/trial?/test/summary.csv  | grep -v overall ; done`

To purge all results:

* `find . -name 'r*.out' | xargs rm -rf`
* `find . -name '.r*.sh*' | xargs rm -f`

To check for errors

* `find . -name '.r*.sh.e' | xargs grep -i error`