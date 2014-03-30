#!/usr/bin/env python

import os
import sys


def handle_dir(dirname):
    with open(os.path.join(dirname, 'Makefile')) as fh:
        in_tss = False
        for ln in fh:
            if ln.startswith('tss:'):
                in_tss = True
            elif in_tss:
                if len(ln.rstrip()) == 0:
                    break
                print >> sys.stderr, '  Found a ts target: %s' % ln.split()[0]


def go():
    for dirname, dirs, files in os.walk('.'):
        if 'Makefile' in files:
            print >> sys.stderr, 'Found a Makefile: %s' % (os.path.join(dirname, 'Makefile'))
            handle_dir(dirname)

go()
