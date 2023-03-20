#!/usr/bin/env python
from __future__ import print_function

"""
python marcc_out.py dry
for dry run: write scripts but doesn't qsub them
python marcc_out.py wet
for normal run: write scripts and also qsub them
"""

import os
import sys
import time
import re


def write_sge(rule, fn, dirname, mem_gb, hours, ncores=8, use_scavenger=False, makefile='Makefile'):
    my_mem_gb, my_hours = mem_gb, hours
    if 'r12' in rule:
        my_mem_gb = int(round(1.5*my_mem_gb))
        my_hours *= 3
        my_hours /= 2
    if 'r12' in rule and 'various_genomes' in rule:
        my_hours *= 3
        my_hours /= 2
    if 'r12' in rule and 'ill_250' in rule:
        my_hours *= 3
        my_hours /= 2
    if 'ill_500' in rule:
        my_hours *= 2
        my_mem_gb = int(round(1.5 * my_mem_gb))
    if '_bwamem' in rule:
        my_mem_gb = max(my_mem_gb, 12)
    if '_snap' in rule:
        my_mem_gb = max(my_mem_gb, 32)
    if 'ts_10m_' in rule:
        my_hours *= 2
    if 'ts_50m_' in rule:
        my_hours *= 4
    pbs_lns = list()
    pbs_lns.append('#!/bin/bash')
    pbs_lns.append('#$ -S /bin/bash')
    pbs_lns.append('#$ -cwd')
    pbs_lns.append('#$ -l mem=%dG' % my_mem_gb)
    pbs_lns.append('#$ -l h_rt=%d:00:00' % my_hours)
    pbs_lns.append('#$ -pe smp %d' % ncores)
    pbs_lns.append('#$ -o ' + fn + '.o')
    pbs_lns.append('#$ -e ' + fn + '.e')
    pbs_lns.append('export QTIP_EXPERIMENTS_HOME=%s' % os.environ['QTIP_EXPERIMENTS_HOME'])
    pbs_lns.append('cd %s' % os.path.abspath(dirname))
    pbs_lns.append('make -f %s %s/DONE' % (makefile, rule))
    with open(os.path.join(dirname, fn), 'w') as ofh:
        ofh.write('\n'.join(pbs_lns) + '\n')

                    
def handle_dir(dirname, re_out, mem_gb, hours, dry_run=True, use_scavenger=False):
    ncores = 8
    with open(os.path.join(dirname, 'Makefile')) as fh:
        in_out = False
        for ln in fh:
            if ln[0] == '#':
                continue
            if ln.startswith('NCORES='):
                ncores = int(ln.split('=')[1])
                assert ncores > 0
            if re_out.match(ln):
                in_out = True
            elif in_out:
                if len(ln.rstrip()) == 0:
                    in_out = False
                else:
                    target = ln.split()[0].split('/')[0]
                    print('  Found a .out target: %s, w/ %d cores' % (target, ncores), file=sys.stderr)
                    target_full = os.path.join(dirname, target)
                    if os.path.exists(os.path.join(target_full, 'DONE')):
                        print('  Skipping target %s because of DONE' % target, file=sys.stderr)
                        continue
                    elif os.path.exists(target_full):
                        # delete it???
                        pass
                    fn = '.' + target + '.sh'
                    write_sge(target, fn, dirname, mem_gb, hours, use_scavenger=use_scavenger, ncores=ncores)
                    print('pushd %s && qsub %s && popd' % (dirname, fn))
                    if not dry_run:
                        os.system('cd %s && qsub %s' % (dirname, fn))
                        time.sleep(0.5)

def write_sge(target, fn, dirname, mem_gb, hours, use_scavenger=False, ncores=8):
    with open(os.path.join(dirname, fn), 'w') as fh:
        fh.write('#!/bin/bash\n')
        fh.write('#$ -S /bin/bash\n')
        fh.write('#$ -cwd\n')
        fh.write('#$ -V\n')
        fh.write('#$ -l h_rt=%d:00:00,h_data=%dG\n' % (hours, mem_gb))
        if use_scavenger:
            fh.write('#$ -l scavenger\n')
        fh.write('#$ -pe smp %d\n' % ncores)
        fh.write('\n')
        fh.write('make %s.out\n' % target)
        fh.write('echo "done" > %s/DONE\n' % target)

def go():
    re_out = re.compile('^outs_[_a-zA-Z01-9]*:.*')
    mem_gb = 8
    hours = 12
    if 'QTIP_EXPERIMENTS_HOME' not in os.environ:
        raise RuntimeError('Must have QTIP_EXPERIMENTS_HOME set')
    for dirname, dirs, files in os.walk('.'):
        if 'Makefile' in files and 'IGNORE' not in files:
            print('Found a Makefile: %s' % (os.path.join(dirname, 'Makefile')), file=sys.stderr)
            handle_dir(dirname, re_out, mem_gb, hours, dry_run=(sys.argv[1] == 'dry' or sys.argv[1] == '--dry'),
                       use_scavenger=len(sys.argv) > 2 and sys.argv[2] == 'scavenger')

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print("pass argument 'dry' for dry run, or 'wet' for normal run")
    else:
        go()

