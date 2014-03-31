#!/usr/bin/env python

__author__ = 'langmead'

import sys

name = None
length = 0
max_len = 10000
first = True
buf = []
for ln in sys.stdin:
    if ln[0] == '>':
        if not first:
            if length < max_len:
                print >> sys.stderr, 'Ref %s has length %d' % (name, length)
            else:
                print '>%s' % name
                print '\n'.join(buf)
                buf = []
        first = False
        length = 0
        name = ln[1:]
    else:
        ln = ln.rstrip()
        length += len(ln)
        buf.append(ln)

if not first:
    if length < max_len:
        print >> sys.stderr, 'Ref %s has length %d' % (name, length)
    else:
        print '>%s' % name
        print '\n'.join(buf)
        buf = []
