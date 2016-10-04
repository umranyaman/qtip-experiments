#!/usr/bin/env python

from __future__ import print_function
from itertools import product


species_list = ['human', 'mouse']
rdlen_list = ['100', '250']
paired_prefix_list = ['r0', 'r1']

headers = ['species', 'rdlen', 'paired', 'type', 'count', 'pct_total', 'pct_error']

print(','.join(headers))
for s, r, p in product(species_list, rdlen_list, paired_prefix_list):
    fn = '_'.join([p, 'mason', s, 'mixture', r]) + '.corstats'
    with open(fn) as fh:
        for ln in fh:
            if ln.startswith('type'):
                continue
            print(','.join([s, r, 'T' if p == 'r1' else 'F', ln.rstrip()]))
