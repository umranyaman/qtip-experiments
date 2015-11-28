#!/usr/bin/env python

"""
Gather together all the results from all the various simulation experiments.
"""

import os
import re
import sys
import glob
import shutil
import logging
from os.path import join


target_re = re.compile('^outs_[_a-zA-Z01-9]*:.*')
fraction = '0.300'
replicate = '1'
sampling_rates = ['0.05', '0.1', '0.2', '0.3', '0.4', '0.5', '1.0']
trials = ['1', '2', '3', '4', '5']


def mkdir_quiet(dr):
    # Create output directory if needed
    import errno
    if not os.path.isdir(dr):
        try:
            os.makedirs(dr)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


def has_done(dr):
    done_fn = join(dr, 'DONE')
    if not os.path.exists(done_fn):
        raise RuntimeError('Directory "%s" does not contain DONE file' % dr)


def copyfiles(fglob, dest, prefix=''):
    assert os.path.isdir(dest) and os.path.exists(dest)
    for fn in glob.glob(fglob):
        shutil.copyfile(fn, join(dest, prefix + os.path.basename(fn)))


def compile_line(fh, combined_target_name, rate, trial, params_fn, summ_fn, first):
    headers = ['name', 'subsampling_rate', 'trial_no']
    values = [combined_target_name, rate, trial]
    for fn in [params_fn, summ_fn]:
        with open(fn, 'r') as fh:
            header = fh.readline().rstrip()
            headers += header.split(',')
            body = fh.readline().rstrip()
            values += body.split(',')
    if first:
        fh.write(','.join(map(str, headers)) + '\n')
    fh.write(','.join(map(str, values)) + '\n')


def handle_dir(dirname, dest_dirname):
    name = os.path.basename(dirname)
    with open(join(dirname, 'Makefile')) as fh:

        in_target = False
        for ln in fh:
            if target_re.match(ln):
                in_target = True
            elif in_target:
                if len(ln.rstrip()) == 0:
                    in_target = False
                else:
                    target = ln.split()[0]
                    target_full = join(dirname, target)
                    has_done(target_full)

                    combined_target_name = name + '_' + target[:-4]
                    odir = join(dest_dirname, combined_target_name)

                    with open(join(odir, 'compiled.tsv')) as fh:

                        first = True
                        for rate in sampling_rates:

                            target_full_s = join(target_full, 'sample' + rate)
                            if not os.path.isdir(target_full_s):
                                raise RuntimeError('Directory "%s" does not exist' % target_full_s)

                            for trial in trials:

                                target_full_st = join(target_full_s, 'trial' + trial)
                                if not os.path.isdir(target_full_st):
                                    raise RuntimeError('Directory "%s" does not exist' % target_full_st)

                                mkdir_quiet(odir)

                                copyfiles(join(target_full_st, 'featimport_*.csv'), odir)
                                params_fn = join(odir, 'params.csv')
                                shutil.copyfile(join(target_full_st, 'params.csv'), params_fn)

                                for tt in ['test', 'training']:

                                    target_full_stt = join(target_full_st, tt)
                                    if not os.path.isdir(target_full_stt):
                                        raise RuntimeError('Directory "%s" does not exist' % target_full_stt)

                                    copyfiles(join(target_full_stt, 'cid*.csv'), odir, tt + '_')
                                    copyfiles(join(target_full_stt, 'cse*.csv'), odir, tt + '_')
                                    copyfiles(join(target_full_stt, 'roc*.csv'), odir, tt + '_')
                                    summ_fn = join(odir, tt + '_summary.csv')
                                    shutil.copyfile(join(target_full_stt, 'summary.csv'), summ_fn)
                                    compile_line(fh, combined_target_name, rate, trial, params_fn, summ_fn, first)
                                    first = False


def go():
    # Set up logger
    logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s',
                        datefmt='%m/%d/%y-%H:%M:%S', level=logging.DEBUG)

    for dirname, dirs, files in os.walk('.'):
        if 'Makefile' in files:
            logging.info('Found a Makefile: %s' % join(dirname, 'Makefile'))
            handle_dir(dirname, 'summary')

    os.system('tar -cvzf summary.tar.gz summary')

go()
