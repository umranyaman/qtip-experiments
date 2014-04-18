#!/usr/bin/env python

import os
import sys
import time


idx = 0
mem_gb = 60


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
    global idx
    with open(os.path.join(dirname, 'Makefile')) as fh:
        in_reads = False
        for ln in fh:
            if ln.startswith('reads:'):
                in_reads = True
            elif in_reads:
                if len(ln.rstrip()) == 0:
                    in_reads = False
                else:
                    target = ln.split()[0]
                    print >> sys.stderr, '  Found a read file: %s' % target
                    target_full = os.path.join(dirname, target)
                    if os.path.exists(target_full):
                        print >> sys.stderr, '  Skipping target %s because target exists' % target
                        continue
                    pbs_lns = list()
                    pbs_lns.append('#PBS -q batch')
                    pbs_lns.append('#PBS -l walltime=60:00')
                    pbs_lns.append('#PBS -j n')
                    for mem_arg in ['pmem', 'vmem', 'pvmem', 'mem']:
                        pbs_lns.append('#PBS -l %s=%dgb' % (mem_arg, mem_gb))
                    pbs_lns.append('export TS_HOME=%s' % os.environ['TS_HOME'])
                    pbs_lns.append('export TS_INDEXES=%s' % os.environ['TS_INDEXES'])
                    pbs_lns.append('export TS_REFS=%s' % os.environ['TS_REFS'])
                    pbs_lns.append('cd %s' % os.path.abspath(dirname))
                    pbs_lns.append('make %s' % target)
                    qsub_dir = '.read_qsubs'
                    mkdir_quiet(qsub_dir)
                    cur_dir = os.getcwd()
                    os.chdir(qsub_dir)
                    qsub_fn = '.%s.%d.sh' % (target, idx)
                    with open(qsub_fn, 'w') as ofh:
                        ofh.write('\n'.join(pbs_lns) + '\n')
                    idx += 1
                    print 'qsub %s' % qsub_fn
                    if not dry_run:
                        os.system('qsub %s' % qsub_fn)
                        time.sleep(0.2)
                    os.chdir(cur_dir)


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
            handle_dir(dirname, dry_run=(len(sys.argv) > 1))

go()
