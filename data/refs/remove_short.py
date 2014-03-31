#!/usr/bin/env python

__author__ = 'langmead'

import sys

name = None
length = 0
max_len = 500
first = True
for ln in sys.stdin:
    if ln[0] == '>':
        if not first and length < max_len:
            print 'Ref %s has length %d' % (name, length)
        first = False
        length = 0
        name = ln[1:].split()[0]
    else:
        length += len(ln.rstrip())
if not first and length < max_len:
    print 'Ref %s has length %d' % (name, length)
