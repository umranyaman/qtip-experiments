#!/usr/bin/env python

"""
Gather together all the results from all the various simulation experiments.
"""

import os
import re
import glob
import logging
from os.path import join


target_re = re.compile('^outs_[_a-zA-Z01-9]*:.*')
fraction = '0.300'
replicate = '1'
sampling_rates = ['0.05', '0.1', '0.2', '0.3', '0.4', '0.5', '1.0']
trials = ['1', '2', '3', '4', '5']


def parse_aligner_local(target):
    # Parsing aligner info
    toks = target.split('_')
    aligner, local = 'bt2', False
    if 'bwamem' in toks[1]:
        aligner = 'bwamem'
        local = True
    if 'snap' in toks[1]:
        aligner = 'snap'
        local = True
    if aligner == 'bt2' and toks[1][-1] == 'l':
        toks[1] = toks[1][:-1]
        local = True
    return aligner, local


def parse_species(target):
    genome = 'hg'
    toks = target.split('_')
    # Parsing reference species info
    if toks[2] == 'mm':
        genome = 'mm'
    if toks[2] == 'zm':
        genome = 'zm'
    return genome


def parse_sim(target):
    sim = 'mason'
    toks = target.split('_')
    # Parsing simulator info:
    if toks[2] == 'wgsim':
        sim = 'wgsim'
    if toks[2] == 'art':
        sim = 'art'
    return sim


def parse_sensitivity(target, aligner):
    toks = target.split('_')
    sensitivity = 's'
    if aligner == 'bt2':
        sensitivity = toks[1][3:]
    return sensitivity


def parse_paired(target):
    return target.startswith('r12')


def parse_readlen(target):
    toks = target.split('_')
    readlen = int(toks[-2])
    if readlen == '50to500':
        readlen = 500
    return readlen


def parse_name_and_target(combined):
    toks = combined.split('_')
    roff = 3
    if toks[2] == 'r0' or toks[2] == 'r12':
        roff = 2
    assert toks[roff] == 'r0' or toks[roff] == 'r12'
    return '_'.join(toks[:roff]), '_'.join(toks[roff:])


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
        os.system('cp -f %s %s' % (fn, join(dest, prefix + os.path.basename(fn))))


def compile_line(ofh, combined_target_name, tt, trial, params_fn, summ_fn, first):
    # TODO: parse name and write some other relevant variables, like whether
    # we're doing local alignment or whether we're aligning pairs
    name, target = parse_name_and_target(combined_target_name)
    aligner, local = parse_aligner_local(target)
    paired = parse_paired(target)
    sim = parse_sim(target)
    readlen = parse_readlen(target)
    sensitivity = parse_sensitivity(target, aligner)
    species = parse_species(target)
    headers = ['name', 'training', 'trial_no', 'aligner', 'local', 'paired',
               'sim', 'readlen', 'sensitivity', 'species']
    values = [name, 'T' if tt == 'training' else 'F', trial, aligner,
              'T' if local else 'F', 'T' if paired else 'F', sim,
              str(readlen), sensitivity, species]
    for fn in [params_fn, summ_fn]:
        with open(fn, 'r') as fh:
            header = fh.readline().rstrip()
            headers += header.split(',')
            body = fh.readline().rstrip()
            values += body.split(',')
    if first:
        ofh.write(','.join(map(str, headers)) + '\n')
    ofh.write(','.join(map(str, values)) + '\n')


def handle_dir(dirname, dest_dirname, ofh, first):
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

                    # Parse some things from the target

                    for rate in sampling_rates:

                        target_full_s = join(target_full, 'sample' + rate)
                        if not os.path.isdir(target_full_s):
                            raise RuntimeError('Directory "%s" does not exist' % target_full_s)

                        for trial in trials:

                            target_full_st = join(target_full_s, 'trial' + trial)
                            if not os.path.isdir(target_full_st):
                                raise RuntimeError('Directory "%s" does not exist' % target_full_st)

                            mkdir_quiet(odir)

                            os.system('cp -f %s %s' % (join(target_full_st, 'featimport_*.csv'), odir))
                            params_fn = join(odir, 'params.csv')
                            os.system('cp -f %s %s' % (join(target_full_st, 'params.csv'), params_fn))

                            for tt in ['test', 'training']:

                                target_full_stt = join(target_full_st, tt)
                                if not os.path.isdir(target_full_stt):
                                    raise RuntimeError('Directory "%s" does not exist' % target_full_stt)

                                copyfiles(join(target_full_stt, 'cid*.csv'), odir, tt + '_')
                                copyfiles(join(target_full_stt, 'cse*.csv'), odir, tt + '_')
                                copyfiles(join(target_full_stt, 'roc*.csv'), odir, tt + '_')
                                summ_fn = join(odir, tt + '_summary.csv')
                                os.system('cp -f %s %s' % (join(target_full_stt, 'summary.csv'), summ_fn))
                                compile_line(ofh, name, combined_target_name, tt, trial, params_fn, summ_fn, first)
                                first = False


def go():
    # Set up logger
    logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s',
                        datefmt='%m/%d/%y-%H:%M:%S', level=logging.DEBUG)

    odir = 'summary'
    first = True
    mkdir_quiet(odir)
    with open(join(odir, 'overall.csv'), 'w') as fh:
        for dirname, dirs, files in os.walk('.'):
            if 'Makefile' in files:
                logging.info('Found a Makefile: %s' % join(dirname, 'Makefile'))
                handle_dir(dirname, odir, fh, first)
                first = False

    os.system('tar -cvzf summary.tar.gz summary')

go()
