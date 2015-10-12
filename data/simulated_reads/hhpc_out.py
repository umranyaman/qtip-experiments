#!/usr/bin/env python

"""
python hhpc_out.py dry

for dry run: write scripts but doesn't qsub them

python hhpc_out.py wet

for normal run: write scripts and also qsub them
"""

import os
import sys
import time
import re


idx = 0
re_out = re.compile('^outs_[_a-zA-Z01-9]*:.*')
mem_gb = 8


def handle_dir(dirname, dry_run=True):
    global idx
    with open(os.path.join(dirname, 'Makefile')) as fh:
        in_out = False
        for ln in fh:
            if ln[0] == '#':
                continue
            if re_out.match(ln):
                in_out = True
            elif in_out:
                if len(ln.rstrip()) == 0:
                    in_out = False
                else:
                    target = ln.split()[0]
                    print >> sys.stderr, '  Found a .out target: %s' % target
                    target_full = os.path.join(dirname, target)
                    if os.path.exists(os.path.join(target_full, 'DONE')):
                        print >> sys.stderr, '  Skipping target %s because of DONE' % target
                        continue
                    elif os.path.exists(target_full):
                        # delete it???
                        pass
                    my_mem_gb = mem_gb
                    if '_bwamem' in target_full:
                        my_mem_gb = 12
                    if '_snap' in target_full:
                        my_mem_gb = 64
                    if '_50M.' in target_full:
                        my_mem_gb = int(my_mem_gb * 1.5)
                    pbs_lns = list()
                    pbs_lns.append('#PBS -q batch')
                    pbs_lns.append('#PBS -l walltime=48:00:00')
                    pbs_lns.append('#PBS -j n')
                    for mem_arg in ['pmem', 'vmem', 'pvmem', 'mem']:
                        pbs_lns.append('#PBS -l %s=%dgb' % (mem_arg, my_mem_gb))
                    pbs_lns.append('export TS_HOME=%s' % os.environ['TS_HOME'])
                    pbs_lns.append('export TS_INDEXES=%s' % os.environ['TS_INDEXES'])
                    pbs_lns.append('export TS_REFS=%s' % os.environ['TS_REFS'])
                    pbs_lns.append('export TMPDIR=/scratch1/langmead-fs1/temp/langmead')
                    pbs_lns.append('cd %s' % os.path.abspath(dirname))
                    pbs_lns.append('if make %s ; then touch %s/DONE ; fi' % (target, target))
                    qsub_basename = '.' + target + '.sh'
                    qsub_fullname = os.path.join(dirname, qsub_basename)
                    with open(qsub_fullname, 'w') as ofh:
                        ofh.write('\n'.join(pbs_lns) + '\n')
                    idx += 1
                    print 'pushd %s && qsub %s && popd' % (dirname, qsub_basename)
                    if not dry_run:
                        os.system('cd %s && qsub %s' % (dirname, qsub_basename))
                        time.sleep(0.5)


def go():
    if 'TS_HOME' not in os.environ:
        raise RuntimeError('Must have TS_HOME set')
    if 'TS_REFS' not in os.environ:
        raise RuntimeError('Must have TS_REFS set')
    if 'TS_INDEXES' not in os.environ:
        raise RuntimeError('Must have TS_INDEXES set')
    for dirname, dirs, files in os.walk('.'):
        if 'Makefile' in files:
            print >> sys.stderr, 'Found a Makefile: %s' % (os.path.join(dirname, 'Makefile'))
            handle_dir(dirname, dry_run=sys.argv[1] == 'dry')

if len(sys.argv) == 1:
    print "pass argument 'dry' for dry run, or 'wet' for normal run"
else:
    go()
