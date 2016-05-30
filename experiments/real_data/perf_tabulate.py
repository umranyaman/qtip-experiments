"""
Creates a table with efficiency metrics (running time, peak memory
footprint) from SLURM output generated from sbatch_align.sh.
"""

from __future__ import print_function

import glob
import sys
from collections import defaultdict

nslurm, nsam = 0, 0
sam_names = defaultdict(int)
tab_wrapped = defaultdict(lambda: defaultdict(int))
to_slurm_wrapped = {}

for fn in glob.glob('slurm-*.out'):
    with open(fn) as fh:
        ln = fh.readline()
        if 'is up to date' in ln:
            continue  # skip trivial job target was already up to date
        if ln.split()[0] != 'python':
            continue

        nsam += 1
        name, t_tandal, t_tandpa, wrappeak, childpeak, t_overall, t_inp = None, 0, 0, 0, 0, 0, 0
        while True:
            ln = fh.readline()
            if len(ln) == 0:
                break
            ln = ln.rstrip()
            if '--vanilla-out' in ln:
                assert name is None
                name = ln.split()[1]
            if 'INFO:Overall' in ln:
                t_overall = float(ln.split()[-1])
            if 'INFO:Aligning input reads' in ln:
                t_inp = float(ln.split()[-1])
            if 'INFO:Aligning tandem reads' in ln and 'paired' not in ln:
                t_tandal = float(ln.split()[-1])
            if 'INFO:Parsing tandem alignments' in ln:
                t_tandpa = float(ln.split()[-1])
            if 'INFO:Peak memory usage (RSS) of Python wrapper' in ln:
                wrappeak = ln.split()[-1]
                assert wrappeak[-2:] == 'GB'
                wrappeak = float(wrappeak[:-2]) * 1024 * 1024 * 1024
            if 'INFO:Peak memory usage (RSS) of children' in ln:
                childpeak = ln.split()[-1]
                assert childpeak[-2:] == 'GB'
                childpeak = float(childpeak[:-2]) * 1024 * 1024 * 1024

        sam_names[name] += 1
        tab_wrapped[name]['wrappeak'] = wrappeak
        tab_wrapped[name]['childpeak'] = childpeak
        tab_wrapped[name]['t_overall'] = t_overall
        tab_wrapped[name]['t_inp'] = t_inp
        tab_wrapped[name]['t_tandal'] = t_tandal
        tab_wrapped[name]['t_tandpa'] = t_tandpa
        to_slurm_wrapped[name] = fn

    nslurm += 1

print('# slurm files: %d' % nslurm, file=sys.stderr)
print('# sam files: %d' % nsam, file=sys.stderr)
for k, v in sorted(sam_names.items()):
    print('  %s: %d' % (k, v), file=sys.stderr)

aln_map = {'bt2': 'Bowtie 2', 'bwa': 'BWA-MEM', 'snap': 'SNAP'}
print('data,aligner,paired,align_time,overall_time,pct_increase_a_to_o,peak_wrapper,peak_children,pct_increase_peak')
for k in sorted(sam_names.keys()):
    wrappeak = tab_wrapped[k]['wrappeak']
    childpeak = tab_wrapped[k]['childpeak']
    wrappct = 0 if childpeak == 0 else (wrappeak * 100.0 / childpeak)
    wrappct = ('+' if wrappct >= 0 else '') + ('%0.3f' % wrappct)
    t_o, t_a = tab_wrapped[k]['t_overall'], tab_wrapped[k]['t_inp']
    t_pct = 0 if t_a == 0 else ((t_o - t_a) * 100.0/t_a)
    t_pct = ('+' if t_pct >= 0 else '') + ('%0.3f' % t_pct)
    srr, aligner, paired, _ = k.split('.')
    srr = srr[:-2]  # chop off trailing _1
    aligner = aln_map[aligner]
    paired = 'T' if paired == 'pair' else 'F'
    print('%s,%s,%s,%.0f,%.0f,%s,%.0f,%.0f,%s' % (srr, aligner, paired, t_a, t_o, t_pct, wrappeak, childpeak, wrappct))
