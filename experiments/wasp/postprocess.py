#!/usr/bin/env python

"""
postprocess.py

Given an original BAM file produced by Qtip, and given a WASP "remap"
FASTQ file, construct a table with one record per remapped read and
with fields giving (a) the original MAPQ, (b) the Qtip-predicted MAPQs
(c) the number of remapped reads that aligned correctly, (d) same, but
incorrectly.

Assumes samtools is in the PATH.

TODO: Only works with Bowtie 2 right now
TODO: Has not been tested with paired-end data yet
"""

from __future__ import print_function
import sys
import os
import re
import subprocess
from collections import defaultdict

join = os.path.join

if 'QTIP_EXPERIMENTS_HOME' not in os.environ:
    raise RuntimeError('QTIP_EXPERIMENTS_HOME must be set')

qtip_exp = os.environ['QTIP_EXPERIMENTS_HOME']

bowtie2_exe = join(qtip_exp, 'software', 'bowtie2', 'bowtie2')


def is_exe(fp):
    return os.path.isfile(fp) and os.access(fp, os.X_OK)


if not is_exe(bowtie2_exe):
    raise RuntimeError('No bowtie2 exe at: "%s"' % bowtie2_exe)


def args_from_bam(bam_fn):
    if not os.path.exists(bam_fn):
        raise RuntimeError('No such BAM as "%s"' % bam_fn)
    proc = subprocess.Popen(['samtools', 'view', '-h', bam_fn], stdout=subprocess.PIPE)
    for ln in proc.stdout:
        if ln[0] != '@':
            raise RuntimeError('Could not parse command line arguments from input bam')
        toks = ln.split('\t')
        if toks[-1].startswith('CL:'):
            cmd = toks[-1][4:-1]
            if cmd.endswith('"'):
                cmd = cmd[:-1]
            cmd = cmd.split(' ')
            cmd = cmd[1:]
            if cmd[0] == '--wrapper':
                cmd = cmd[2:]
            return cmd


def remove_args(bt2_args, exclude):
    for i in range(len(bt2_args)-1, -1, -1):
        if bt2_args[i] in exclude:
            bt2_args = bt2_args[:i] + bt2_args[i+2:]
    return bt2_args


def align_fastq(fastq_fn, bt2_args, threads, ofn):
    if not os.path.exists(fastq_fn):
        raise RuntimeError('No such FASTQ as "%s"' % fastq_fn)
    ofn += '.sam'
    ofn = join(os.path.dirname(ofn), '.' + os.path.basename(ofn))
    cmd = [bowtie2_exe] + remove_args(bt2_args, ['-U', '-1', '-2', '-S', '-p'])
    cmd += ['-U', fastq_fn, '-S', ofn, '-p', str(threads)]
    cmd = ' '.join(cmd)
    print('  Bowtie2 command: ' + cmd, file=sys.stderr)
    ret = os.system(cmd)
    if ret != 0:
        raise RuntimeError('Bowite2 process returned %d' % ret)
    return ofn


# Scan the resulting BAM together with the original BAM
def scan_remapped_bam(remapped_sam_fn):
    hist = defaultdict(lambda: [0, 0])  # read name -> [# correct, # incorrect]
    with open(remapped_sam_fn) as fh:
        for ln in fh:
            if ln[0] == '@':
                continue
            toks = ln.split('\t')
            words = toks[0].split(".")

            # ********************************************************
            # The following code is borrowed and adapted from WASP.
            # WASP is by Bryce van de Geijn, Graham McVicker, Yoav
            # Gilad, & Jonathan Pritchard and is available here:
            # https://github.com/bmvdgeijn/WASP
            # ********************************************************

            if len(words) < 4:
                raise ValueError("expected read names to be formatted "
                                 "like <orig_name>.<coordinate>."
                                 "<read_number>.<total_read_number> but got "
                                 "%s" % toks[0])

            # token separator '.' can potentially occur in
            # original read name, so if more than 4 tokens,
            # assume first tokens make up original read name
            coord_str, num_str, total_str = words[-3:]
            orig_name = ".".join(words[0:-3])
            flags, pos = int(toks[1]), int(toks[3])
            next_reference_start = int(toks[7])
            correct_map = int(coord_str) == pos

            if '-' in coord_str:
                # paired end read, coordinate gives expected positions for each end
                c1, c2 = coord_str.split("-")

                if (flags & 3) != 3:
                    continue  # not paired or not proper pair

                pos1, pos2 = int(c1), int(c2)

                # only use left end of reads, but check that right end is in
                # correct location
                if pos < next_reference_start:
                    correct_map = pos1 == pos + 1 and pos2 == next_reference_start + 1
                else:
                    continue  # this is right end of read

            hist[orig_name][0 if correct_map else 1] += 1

    return hist


# Output a file with one line per input alignment:
# - Original MAPQ
# - Predicted MAPQ
# - # derived reads that aligned correctly
# - # derived reads that aligned incorrectly
def tabulate(bam_fn, out_fn, hist):
    mapq_re = re.compile('Zm:[iZ]:([0-9]+)')
    with open(out_fn, 'wb') as ofh:
        proc = subprocess.Popen(['samtools', 'view', '-h', bam_fn], stdout=subprocess.PIPE)
        for ln in proc.stdout:
            if ln[0] == '@':
                continue
            toks = ln.split('\t')
            qname = toks[0]
            if qname in hist:
                mapq = toks[4]
                cor, incor = hist[qname]
                orig_mapq = mapq_re.search(ln)
                if orig_mapq is None:
                    raise RuntimeError('Could not parse original mapq from this line: ' + ln)
                orig_mapq = orig_mapq.group(1)
                print(','.join([mapq, orig_mapq, str(cor), str(incor)]), file=ofh)


def add_args(parser):
    parser.add_argument('--bam', metavar='path', type=str, required=True, help='Input BAM file')
    parser.add_argument('--fastq', metavar='path', type=str, required=True, help='WASP "remap" FASTQ file')
    parser.add_argument('--output', metavar='path', type=str, required=True, help='Place output table here')
    parser.add_argument('--threads', type=int, default=16, help='Use -p with Bowtie 2')


def go(args):
    print('Getting arguments from BAM', file=sys.stderr)
    bt2_args = args_from_bam(args.bam)
    print('  Bowtie2 arguments: ' + str(bt2_args), file=sys.stderr)
    print('Aligning overlap FASTQ', file=sys.stderr)
    remap_sam_fn = align_fastq(args.fastq, bt2_args, args.threads, args.output)
    print('Scanning alignments', file=sys.stderr)
    hist = scan_remapped_bam(remap_sam_fn)
    print('Tabulating', file=sys.stderr)
    tabulate(args.bam, args.output, hist)


def go_profile(args):
    if args.profile:
        import cProfile
        cProfile.run('go(args)')
    else:
        go(args)


if __name__ == "__main__":

    import argparse

    _parser = argparse.ArgumentParser(
        description='Given an original BAM file produced by Qtip, and '
                    'given a WASP remap FASTQ file, construct a table '
                    'giving original and predicted MAPQs together with '
                    'the number of remapped reads that aligned '
                    'correctly and uncorrectly')

    add_args(_parser)
    _parser.add_argument('--profile', action='store_const', const=True, default=False, help='Print profiling info')
    _args = _parser.parse_args(sys.argv[1:])

    go_profile(_args)
