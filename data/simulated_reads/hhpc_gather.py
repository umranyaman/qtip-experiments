#!/usr/bin/env python

import os
import re
import shutil
import logging


pred_re = re.compile('preds[_a-zA-Z01-9]*:.*')
fraction = '0.300'
replicate = '1'


def mkdir_quiet(dr):
    # Create output directory if needed
    import errno
    if not os.path.isdir(dr):
        try:
            os.makedirs(dr)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


def handle_dir(dirname):
    with open(os.path.join(dirname, 'Makefile')) as fh:
        in_reads = False
        for ln in fh:
            if pred_re.match(ln):
                in_reads = True
            elif in_reads:
                if len(ln.rstrip()) == 0:
                    in_reads = False
                else:
                    target = ln.split()[0]
                    name = os.path.basename(dirname)
                    target_full = os.path.join(dirname, target)

                    # subsampling png
                    subsampling_png_src_fn = 'subsampling_series.png'
                    subsampling_png_dst_fn = '%s_%s_subsampling_series.png' % (name, target)
                    subsampling_png_full = os.path.join(target_full, subsampling_png_src_fn)
                    subsampling_tsv_out_dir = os.path.join('summary', 'subsampling_plots')
                    mkdir_quiet(subsampling_tsv_out_dir)
                    subsampling_png_out = os.path.join(subsampling_tsv_out_dir, subsampling_png_dst_fn)
                    if os.path.exists(subsampling_png_full):
                        shutil.copyfile(subsampling_png_full, subsampling_png_out)
                    else:
                        logging.warning('Could not find source file "%s"' % subsampling_png_full)

                    # subsampling tsvs
                    subsampling_tsv_out_dir = os.path.join('summary', 'subsampling_tables')
                    mkdir_quiet(subsampling_tsv_out_dir)
                    for i in xrange(1, 6):
                        subsampling_tsv_src_fn = 'subsampling_series_%d.tsv' % i
                        subsampling_tsv_dst_fn = '%s_%s_subsampling_series_%d.tsv' % (name, target, i)
                        subsampling_tsv_full = os.path.join(target_full, subsampling_tsv_src_fn)
                        subsampling_tsv_out = os.path.join(subsampling_tsv_out_dir, subsampling_tsv_dst_fn)
                        if os.path.exists(subsampling_tsv_full):
                            shutil.copyfile(subsampling_tsv_full, subsampling_tsv_out)
                        else:
                            logging.warning('Could not find source file "%s"' % subsampling_tsv_full)

                    # ROC tables
                    roc_src_fn = 'roc_table.tsv'
                    roc_src_orig_fn = 'roc_table_orig.tsv'
                    roc_tsv_full = os.path.join(target_full, 'subsampled', fraction, replicate, roc_src_fn)
                    roc_orig_tsv_full = os.path.join(target_full, 'subsampled', fraction, replicate, roc_src_orig_fn)
                    roc_tsv_out_dir = os.path.join('summary', 'roc_table')
                    mkdir_quiet(roc_tsv_out_dir)
                    roc_dst_fn = '%s_%s_roc_table.tsv' % (name, target)
                    roc_tsv_out = os.path.join(roc_tsv_out_dir, roc_dst_fn)
                    roc_dst_orig_fn = '%s_%s_roc_table_orig.tsv' % (name, target)
                    roc_tsv_orig_out = os.path.join(roc_tsv_out_dir, roc_dst_orig_fn)
                    for fn, ofn in zip([roc_tsv_full, roc_orig_tsv_full],
                                       [roc_tsv_out, roc_tsv_orig_out]):
                        if os.path.exists(fn):
                            shutil.copyfile(fn, ofn)
                        else:
                            logging.warning('Could not find source file "%s"' % fn)


def go():
    # Set up logger
    logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s',
                        datefmt='%m/%d/%y-%H:%M:%S', level=logging.DEBUG)

    for dirname, dirs, files in os.walk('.'):
        if 'Makefile' in files:
            logging.info('Found a Makefile: %s' % os.path.join(dirname, 'Makefile'))
            handle_dir(dirname)

    os.system('tar -cvzf summary.tar.gz summary')

go()
