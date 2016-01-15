from __future__ import print_function

import glob
from collections import defaultdict

nslurm, nsam, nplain = 0, 0, 0
sam_names = defaultdict(int)
plain_names = defaultdict(int)
tab_plain = defaultdict(lambda: defaultdict(int))
tab_wrapped = defaultdict(lambda: defaultdict(int))

for fn in glob.glob('slurm-*.out'):
    with open(fn) as fh:
        ln = fh.readline()
        if ln.split()[0] == 'python':
            nsam += 1
            name, maxmss = None, 0
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

            sam_names[name] += 1
            #assert maxmss > 0
            tab_wrapped[name]['maxmss'] = maxmss

        else:
            nplain += 1
            name, maxmss = None, 0
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

            plain_names[name] += 1
            #assert maxmss > 0
            tab_plain[name]['maxmss'] = maxmss

    nslurm += 1

print('# slurm files: %d' % nslurm)
print('# sam files: %d' % nsam)
for k, v in sorted(sam_names.items()):
    print('  %s: %d' % (k, v))
print('# plain files: %d' % nplain)
for k, v in sorted(plain_names.items()):
    print('  %s: %d' % (k, v))

for k in sorted(sam_names.keys()):
    maxmss = tab_wrapped[k]['maxmss']
    maxmss_plain = tab_plain[k]['maxmss']
    diff = maxmss - maxmss_plain
    diff_pct = 0 if maxmss_plain == 0 else (100.0 * float(diff) / maxmss_plain)
    diff_pct = ('+' if diff_pct >= 0 else '') + ('%0.2f' % diff_pct)
    print('%s: %d %d %s' % (k, maxmss_plain, maxmss, diff_pct))
