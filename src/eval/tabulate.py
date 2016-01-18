from __future__ import print_function

import glob
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

print('# slurm files: %d' % nslurm)
print('# sam files: %d' % nsam)
for k, v in sorted(sam_names.items()):
    print('  %s: %d' % (k, v))

for k in sorted(sam_names.keys()):
    wrappeak = tab_wrapped[name]['wrappeak']
    childpeak = tab_wrapped[name]['childpeak']
    wrappct = 0 if childpeak == 0 else (wrappeak * 100.0 / childpeak)
    wrappct = ('+' if wrappct >= 0 else '') + ('%0.3f' % wrappct)
    t_o, t_a = tab_wrapped[name]['t_overall'], tab_wrapped[name]['t_inp']
    t_pct = 0 if t_a == 0 else ((t_o - t_a) * 100.0/t_a)
    t_pct = ('+' if t_pct >= 0 else '') + ('%0.3f' % t_pct)
    print('%s: %.0f,%.0f,%s %.0f,%.0f,%s %s' % (k, t_a, t_o, t_pct, wrappeak, childpeak, wrappct, to_slurm_wrapped[k]))
