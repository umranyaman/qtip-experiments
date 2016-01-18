from __future__ import print_function

import glob
from collections import defaultdict

nslurm, nsam, nplain = 0, 0, 0
sam_names = defaultdict(int)
plain_names = defaultdict(int)
tab_plain = defaultdict(lambda: defaultdict(int))
tab_wrapped = defaultdict(lambda: defaultdict(int))
to_slurm_wrapped = {}
to_slurm_plain = {}

for fn in glob.glob('slurm-*.out'):
    with open(fn) as fh:
        ln = fh.readline()
        if 'is up to date' in ln:
            continue  # skip trivial job target was already up to date
        if ln.split()[0] == 'python':
            nsam += 1
            name, maxmss, wallcl, tandal, tandpa, wrappeak, childpeak = None, 0, 0, 0, 0, 0, 0
            while True:
                ln = fh.readline()
                if len(ln) == 0:
                    break
                ln = ln.rstrip()
                if '--vanilla-out' in ln:
                    assert name is None
                    name = ln.split()[1]
                if 'Maximum resident set size' in ln:
                    maxmss = int(ln.split()[-1])
                if 'Elapsed (wall clock) time' in ln:
                    toks = map(int, ln.split()[-1].split(':'))
                    wallcl = toks[0]*60*60 + toks[1]*60 + toks[2]
                if 'INFO:Aligning tandem reads' in ln and 'paired' not in ln:
                    tandal = float(ln.split()[-1])
                if 'INFO:Parsing tandem alignments' in ln:
                    tandpa = float(ln.split()[-1])
                if 'INFO:Peak memory usage (RSS) of Python wrapper' in ln:
                    wrappeak = ln.split()[-1]
                    assert wrappeak[-2:] == 'GB'
                    wrappeak = float(wrappeak[:-2]) * 1024 * 1024 * 1024
                if 'INFO:Peak memory usage (RSS) of children' in ln:
                    childpeak = ln.split()[-1]
                    assert childpeak[-2:] == 'GB'
                    childpeak = float(childpeak[:-2]) * 1024 * 1024 * 1024

            sam_names[name] += 1
            tab_wrapped[name]['maxmss'] = maxmss
            tab_wrapped[name]['wrappeak'] = wrappeak
            tab_wrapped[name]['childpeak'] = childpeak
            tab_wrapped[name]['wallcl'] = wallcl
            tab_wrapped[name]['tandal'] = tandal
            tab_wrapped[name]['tandpa'] = tandpa
            to_slurm_wrapped[name] = fn

        else:
            nplain += 1
            name, maxmss, wallcl = None, 0, 0
            while True:
                ln = fh.readline()
                if len(ln) == 0:
                    break
                ln = ln.rstrip()
                if ln.endswith('.plain'):
                    tok = ln.split()[-1][:-6]
                    assert name is None or name == tok
                    name = tok
                if 'Maximum resident set size' in ln:
                    maxmss = int(ln.split()[-1])
                if 'Elapsed (wall clock) time' in ln:
                    toks = map(int, ln.split()[-1].split(':'))
                    wallcl = toks[0]*60*60 + toks[1]*60 + toks[2]

            plain_names[name] += 1
            tab_plain[name]['maxmss'] = maxmss
            tab_plain[name]['wallcl'] = wallcl
            to_slurm_plain[name] = fn

    nslurm += 1

print('# slurm files: %d' % nslurm)
print('# sam files: %d' % nsam)
for k, v in sorted(sam_names.items()):
    print('  %s: %d' % (k, v))
print('# plain files: %d' % nplain)
for k, v in sorted(plain_names.items()):
    print('  %s: %d' % (k, v))

for k in sorted(sam_names.keys()):
    wallcl = tab_wrapped[k]['wallcl']
    wallcl_plain = tab_plain[k]['wallcl']
    diffcl = wallcl - wallcl_plain
    diffcl_pct = 0 if wallcl_plain == 0 else (100.0 * float(diffcl) / wallcl_plain)
    diffcl_pct = ('+' if diffcl_pct >= 0 else '') + ('%0.2f' % diffcl_pct)
    difftandem = tab_wrapped[k]['tandal'] + tab_wrapped[k]['tandpa']
    difftandem_pct = 0 if wallcl_plain == 0 else (100.0 * float(difftandem) / wallcl_plain)
    difftandem_pct = ('+' if difftandem_pct >= 0 else '') + ('%0.2f' % difftandem_pct)
    wrappeak = tab_wrapped[name]['wrappeak']
    childpeak = tab_wrapped[name]['childpeak']
    wrappct = 0 if childpeak == 0 else (wrappeak * 100.0 / childpeak)
    wrappct = ('+' if wrappct >= 0 else '') + ('%0.2f' % wrappct)
    print('%s: %d,%d,%s,%s %.0f,%.0f,%s %s %s' % (k, wallcl, wallcl_plain, diffcl_pct, difftandem_pct,
                                                     wrappeak, childpeak, wrappct,
                                                     to_slurm_plain[k], to_slurm_wrapped[k]))
