#!/usr/bin/env python

"""
Output should have columns for:

- aligner
- local/global
- strict/loose
- unpaired/paired
- roc_round
- roc_orig

Once the ROCs are parsed, R code can recreate RCA & RCE
"""

from __future__ import print_function
import glob


def roc_file_to_string(roc_fn, inner_sep=':', outer_sep=';'):
    """ Convert a file with a ROC table into a string with one line per ROC row """
    fields = []
    with open(roc_fn) as fh:
        for ln in fh:
            mapq, cor, incor, _, _ = ln.rstrip().split(',')
            if mapq == "mapq":
                continue
            fields.append(inner_sep.join([mapq, cor, incor]))
    return outer_sep.join(fields)


def parse_fn(fn):
    """ ERR050083_1.bwamem-qtip_loose_dec.unp.roc.csv """
    ptoks = fn.split('.')
    data = ptoks[0][:-2]
    paired = 'T' if ptoks[2] == 'pair' else 'F'
    aligner = ptoks[1].split('-')[0]
    local = 'T' if aligner != 'bt2' else 'F'
    loose = 'T' if ptoks[1].split('_')[1] == 'loose' else 'F'
    return data, paired, aligner, local, loose


print(','.join(["data", "paired", "aligner", "local", "loose", "roc_round", "roc_orig"]))
for fn in glob.glob('*-qtip_*int.*.roc.csv'):
    roc_str = roc_file_to_string(fn)
    roc_orig_str = roc_file_to_string(fn.replace('_int.', '_orig.'))
    data, paired, aligner, local, loose = parse_fn(fn)
    print(','.join([data, paired, aligner, local, loose, roc_str, roc_orig_str]))
