#!/usr/bin/env python

"""
Gathers results from real_data "correctish" experiments.  Descends into
experimental subdirectories and parses the Makefiles it finds.

Outputs:
 - "overall.csv" with one big table of summary measures
 - "summary" subdirectory with lots of raw results compiled into a directory
   structure
   + Subdirectories correspond to simulation experiments and contain:
     - for each sampling rate:
       + for each trial:
         - featimport_*.csv -- feature importances for each alignment type
         - params.csv -- feature importances for each model
         - Subdirectories for training/test, each with:
           + roc.csv -- ROC table
           + summary.csv -- summarizes data, model fit
"""

from __future__ import print_function
import sys
import os
import logging


def roc_file_to_string(roc_fn, inner_sep=':', outer_sep=';'):
    """ Convert a file with a ROC table into a string with one line per ROC row """
    fields = []
    with open(roc_fn) as fh:
        for ln in fh:
            cor, _, _, incor, mapq, _, _, _, _, _ = ln.rstrip().split(',')
            if mapq == "mapq":
                continue
            fields.append(inner_sep.join([mapq, cor, incor]))
    return outer_sep.join(fields)


# new_ERR050083_1.bt2.unp.sam
def compile_line(ofh, name, trial, al, training, params_fn, summ_fn, roc_round_fn, roc_orig_fn, first):
    headers = ['name', 'trial_no', 'aligner', 'training']
    values = [name, trial, al, 'T' if training else 'F']
    for fn in [params_fn, summ_fn]:
        with open(fn, 'r') as fh:
            header = fh.readline().rstrip()
            headers += header.split(',')
            body = fh.readline().rstrip()
            values += body.split(',')
    # Add ROCs; these are big long strings
    headers.extend(['roc_round', 'roc_orig'])
    values.extend([roc_file_to_string(roc_round_fn),
                   roc_file_to_string(roc_orig_fn)])
    if first:
        ofh.write(','.join(map(str, headers)) + '\n')
    ofh.write(','.join(map(str, values)) + '\n')


def go():

    # Set up logger
    logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s',
                        datefmt='%m/%d/%y-%H:%M:%S', level=logging.DEBUG)

    out_fn = 'overall.csv'
    first = True
    with open(out_fn, 'w') as ofh:
        for al in ['bt2', 'bwa', 'snap']:
            for trial in '0123456789':
                for samp in ['ERR050082', 'ERR050083']:
                    trial_dr = os.path.join('new_%s_1.%s.unp.sam' % (samp, al), 'trial%d' % trial)
                    params_fn = os.path.join(trial_dr, 'params.csv')
                    trial_tr_dr = os.path.join(trial_dr, 'train')
                    trial_te_dr = os.path.join(trial_dr, 'test')
                    tr_summary = os.path.join(trial_tr_dr, 'summary.csv')
                    te_summary = os.path.join(trial_te_dr, 'summary.csv')
                    tr_roc = os.path.join(trial_tr_dr, 'roc_round.csv')
                    te_roc = os.path.join(trial_te_dr, 'roc_round.csv')
                    tr_roc_orig = os.path.join(trial_tr_dr, 'roc_orig.csv')
                    te_roc_orig = os.path.join(trial_te_dr, 'roc_orig.csv')
                    compile_line(ofh, samp, trial, al, True,  params_fn, tr_summary, tr_roc, tr_roc_orig, first)
                    first = False
                    compile_line(ofh, samp, trial, al, False, params_fn, te_summary, te_roc, te_roc_orig, first)


if '--slurm' in sys.argv:
    script_fn = '.gather.sh'
    gather_args = ''
    if '--experiment' in sys.argv:
        exp_name = sys.argv[sys.argv.index('--experiment')+1]
        script_fn = '.gather_%s.sh' % exp_name
        gather_args = '--experiment ' + exp_name
    my_hours = 4
    with open(script_fn, 'w') as ofh:
        print("#!/bin/bash -l", file=ofh)
        print("#SBATCH", file=ofh)
        print("#SBATCH --nodes=1", file=ofh)
        print("#SBATCH --mem=4G", file=ofh)
        if '--scavenger' in sys.argv:
            print('#SBATCH --partition=scavenger', file=ofh)
            print('#SBATCH --qos=scavenger', file=ofh)
        else:
            print('#SBATCH --partition=shared', file=ofh)
        print('#SBATCH --time=%d:00:00' % my_hours, file=ofh)
        print('#SBATCH --output=%s.o' % script_fn, file=ofh)
        print('#SBATCH --error=%s.e' % script_fn, file=ofh)
        print('python gather.py %s' % gather_args, file=ofh)
    print('sbatch ' + script_fn)
else:
    go()
