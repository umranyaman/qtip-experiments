#!/usr/bin/env python
from __future__ import print_function

"""
python marcc_reads.py dry
for dry run: write scripts but doesn't qsub them
python marcc_reads.py wet
for normal run: write scripts and also qsub them
"""

import os
import sys
import time


idx = 0
mem_gb = 64
hours = 16
jobs = 0


def mkdir_quiet(dr):
    # Create output directory if needed
    import errno
    if not os.path.isdir(dr):
        try:
            os.makedirs(dr)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


def handle_dir(dirname, dry_run=True):
    global idx, jobs
    with open(os.path.join(dirname, 'Makefile')) as fh:
        in_reads = False
        for ln in fh:
            if ln[0] == '#':
                continue
            if ln.startswith('reads:'):
                in_reads = True
            elif in_reads:
                if len(ln.rstrip()) == 0:
                    in_reads = False
                else:
                    target = ln.split()[0]
                    print('  Found a read target: %s' % target, file=sys.stderr)
                    target_full = os.path.join(dirname, target)
                    if os.path.exists(target_full):
                        print('  Skipping target %s because target exists' % target, file=sys.stderr)
                        continue
                    my_mem_gb, my_hours = mem_gb, hours
                    qsub_basename = '.' + target + '.sh'
                    qsub_lns = list()
                    qsub_lns.append('#!/bin/bash')
                    qsub_lns.append('#$ -l h_rt=%d:00:00' % my_hours)
                    qsub_lns.append('#$ -l h_vmem=%dG' % my_mem_gb)
                    qsub_lns.append('#$ -o ' + qsub_basename + '.o')
                    qsub_lns.append('#$ -e ' + qsub_basename + '.e')
                    qsub_lns.append('export QTIP_EXPERIMENTS_HOME=%s' % os.environ['QTIP_EXPERIMENTS_HOME'])
                    qsub_lns.append('cd %s' % os.path.abspath(dirname))
                    qsub_lns.append('make %s' % target)
                    qsub_fullname = os.path.join(dirname, qsub_basename)
                    with open(qsub_fullname, 'w') as ofh:
                        ofh.write('\n'.join(qsub_lns) + '\n')
                    idx += 1
                    print('pushd %s && qsub %s && popd' % (dirname, qsub_basename))
                    jobs += 1
                    if not dry_run:
                        os.system('cd %s && qsub %s' % (dirname, qsub_basename))
                        time.sleep(0.5)


def go():
    if 'QTIP_EXPERIMENTS_HOME' not in os.environ:
        raise RuntimeError('Must have QTIP_EXPERIMENTS_HOME set')
    for dirname, dirs, files in os.walk('.'):
        if 'Makefile' in files and 'IGNORE' not in files:
            print('Found a Makefile: %s' % (os.path.join(dirname, 'Makefile')), file=sys.stderr)
            handle_dir(dirname, dry_run=sys.argv[1] == 'dry')
    print('Composed %d jobs' % jobs, file=sys.stderr)

if len(sys.argv) == 1:
    print("pass argument 'dry' for dry run, or 'wet' for normal run")
else:
    go()
