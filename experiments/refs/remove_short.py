#!/usr/bin/env python

from __future__ import print_function
import sys

name = None
length = 0
max_len = 10000
first = True
buf = []
short_fh = None
if len(sys.argv) > 1:
    short_fh = open(sys.argv[1], 'w')


def handle_entry():
    filt = length < max_len
    print('length(%s) = %d, %s' % (name, length, 'TOO SHORT' if filt else 'pass'), file=sys.stderr)
    if filt:
        if short_fh is not None:
            print('>' + name, file=short_fh)
            print('\n'.join(buf), file=short_fh)
    else:
        print('>' + name)
        print('\n'.join(buf))


for ln in sys.stdin:
    if ln[0] == '>':
        if not first:
            handle_entry()
            buf = []
        first = False
        length = 0
        name = ln[1:].rstrip()
    else:
        ln = ln.rstrip()
        length += len(ln)
        buf.append(ln)

if not first:
    handle_entry()

if short_fh is not None:
    short_fh.close()
