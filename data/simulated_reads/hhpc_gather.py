#!/usr/bin/env python

import os
import sys
import re
import shutil


pred_re = re.compile('preds[_a-zA-Z01-9]*:.*')


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
                    subsampling_png_out_dir = os.path.join('summary', 'subsampling_plots')
                    mkdir_quiet(subsampling_png_out_dir)
                    subsampling_png_out = os.path.join(subsampling_png_out_dir, subsampling_png_dst_fn)
                    shutil.copyfile(subsampling_png_full, subsampling_png_out)


def go():
    for dirname, dirs, files in os.walk('.'):
        if 'Makefile' in files:
            print >> sys.stderr, 'Found a Makefile: %s' % (os.path.join(dirname, 'Makefile'))
            handle_dir(dirname)

go()
