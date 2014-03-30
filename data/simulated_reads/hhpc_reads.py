#!/usr/bin/env python

import os
import sys


idx = 0


def handle_dir(dirname):
    global idx
    with open(os.path.join(dirname, 'Makefile')) as fh:
        in_reads = False
        for ln in fh:
            if ln.startswith('reads:'):
                in_reads = True
            elif in_reads:
                if len(ln.rstrip()) == 0:
                    break
                target = ln.split()[0]
                print >> sys.stderr, '  Found a read file: %s' % target
                pbs_lns = list()
                pbs_lns.append('#PBS -q batch')
                pbs_lns.append('#PBS -l walltime=15:00')
                pbs_lns.append('#PBS -j n')
                pbs_lns.append('#PBS -l pmem=12gb')
                pbs_lns.append('cd %s' % dirname)
                pbs_lns.append('make %s' % target)
                qsub_fn = '.%s.%d.sh' % (target, idx)
                with open(qsub_fn, 'w') as ofh:
                    ofh.write('\n'.join(pbs_lns) + '\n')
                idx += 1
                print 'qsub %s' % qsub_fn


def go():
    for dirname, dirs, files in os.walk('.'):
        if 'Makefile' in files:
            print >> sys.stderr, 'Found a Makefile: %s' % (os.path.join(dirname, 'Makefile'))
            handle_dir(dirname)

go()
