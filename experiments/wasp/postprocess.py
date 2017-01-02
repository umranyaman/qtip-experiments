#!/usr/bin/env python

"""
postprocess.py

Given an original BAM file produced by Qtip, and given a WASP "remap"
FASTQ file, construct a table with one record per remapped read and
with fields giving (a) the original MAPQ, (b) the Qtip-predicted MAPQs
(c) the number of remapped reads that aligned correctly, (d) same, but
incorrectly.

Assumes samtools is in the PATH.
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
bwa_exe = join(qtip_exp, 'software', 'bwa', 'bwa')
snap_exe = join(qtip_exp, 'software', 'snap', 'snap-aligner')


def is_exe(fp):
    return os.path.isfile(fp) and os.access(fp, os.X_OK)


def args_from_bam(bam_fn):
    if not os.path.exists(bam_fn):
        raise RuntimeError('No such BAM as "%s"' % bam_fn)
    proc = subprocess.Popen(['samtools', 'view', '-h', bam_fn], stdout=subprocess.PIPE)
    for ln in proc.stdout:
        if ln[0] != '@':
            raise RuntimeError('Could not parse command line arguments from input bam')
        cmd, myid = None, None
        for tok in ln.rstrip().split('\t'):
            if tok.startswith('ID:'):
                myid = tok[3:]
            if tok.startswith('CL:'):
                cmd = tok[3:]
                if cmd.startswith('"'):
                    cmd = cmd[1:]
                if cmd.endswith('"'):
                    cmd = cmd[:-1]
                cmd = cmd.split(' ')
                # snap doesn't put the binary in the CL: field
                if myid != 'SNAP':
                    cmd = cmd[1:]
                if cmd[0] == '--wrapper':
                    cmd = cmd[2:]
            if cmd is not None and myid is not None:
                proc.terminate()
                return myid, cmd
    raise RuntimeError('should not be here')


def remove_args(bt2_args, exclude, num_to_remove=1):
    for i in range(len(bt2_args)-1, -1, -1):
        if bt2_args[i] in exclude:
            bt2_args = bt2_args[:i] + bt2_args[i+1+num_to_remove:]
    return bt2_args


def align_fastq_bwa(fastq1_fn, fastq2_fn, bwa_args, threads, ofn):
    while bwa_args[-1].endswith('.fastq') or bwa_args[-1].endswith('.fq') or \
            bwa_args[-1].endswith('.fastq.gz') or bwa_args[-1].endswith('.fq.gz'):
        bwa_args = bwa_args[:-1]
    cmd = [bwa_exe] + remove_args(bwa_args, ['-t']) + ['-t', str(threads), fastq1_fn]
    if fastq2_fn is not None:
        cmd.append(fastq2_fn)
    cmd.extend(['>', ofn])
    cmd = ' '.join(cmd)
    print('  BWA command: ' + cmd, file=sys.stderr)
    ret = os.system(cmd)
    if ret != 0:
        raise RuntimeError('BWA process returned %d' % ret)
    return ofn


def align_fastq_bowtie2(fastq1_fn, fastq2_fn, bt2_args, threads, ofn):
    cmd = [bowtie2_exe] + remove_args(bt2_args, ['-U', '-1', '-2', '-S', '-p'])
    if fastq2_fn is None:
        cmd.extend(['-U', fastq1_fn])
    else:
        cmd.extend(['-1', fastq1_fn, '-2', fastq2_fn])
    cmd.extend(['-S', ofn, '-p', str(threads)])
    cmd = ' '.join(cmd)
    print('  Bowtie2 command: ' + cmd, file=sys.stderr)
    ret = os.system(cmd)
    if ret != 0:
        raise RuntimeError('Bowite2 process returned %d' % ret)
    return ofn


def align_fastq_snap(fastq1_fn, fastq2_fn, snap_args, threads, ofn):
    assert '-o' in snap_args
    cmd = [snap_exe] + remove_args(remove_args(snap_args, ['-sam', '-t']), ['-fastq', '-compressedFastq'], 1 if fastq2_fn is None else 2)
    assert '-o' in snap_args
    format_arg = '-compressedFastq' if fastq1_fn.endswith('.gz') else '-fastq'
    if fastq2_fn is None:
        cmd = cmd[0:3] + [format_arg, fastq1_fn] + cmd[3:]
    else:
        cmd = cmd[0:3] + [format_arg, fastq1_fn, fastq2_fn] + cmd[3:]
    oidx = cmd.index('-o')
    cmd = cmd[0:oidx+1] + ['-sam', ofn] + cmd[oidx+1:]
    cmd.extend(['-t', str(threads)])
    cmd.extend(['-xf', '4.0'])
    cmd = ' '.join(cmd)
    print('  SNAP command: ' + cmd, file=sys.stderr)
    ret = os.system(cmd)
    if ret != 0:
        raise RuntimeError('SNAP process returned %d' % ret)
    return ofn


def print_update(nsecondary, nnot_proper_pair, ncor, nincor, nleft, nright):
    print('  ---', file=sys.stderr)
    print('  # secondary: %d' % nsecondary, file=sys.stderr)
    print('  # not paired or proper: %d' % nnot_proper_pair, file=sys.stderr)
    print('  # left: %d, # right: %d' % (nleft, nright), file=sys.stderr)
    if ncor + nincor > 0:
        pct = float(ncor)*100.0/(ncor+nincor)
        print('  # correct: %d (%0.2f%%), # incorrect: %d' % (ncor, pct, nincor), file=sys.stderr)


# Scan the resulting BAM together with the original BAM
def scan_remapped_bam(remapped_sam_fn, keep=False):
    hist = defaultdict(lambda: [0, 0])  # read name -> [# correct, # incorrect]
    n, nsecondary, nnot_proper_pair, nleft, nright = 0, 0, 0, 0, 0
    ncor, nincor = 0, 0
    nival = 10
    nival_fact = 1.25
    with open(remapped_sam_fn) as fh:
        for ln in fh:
            n += 1
            if n == nival:
                print_update(nsecondary, nnot_proper_pair, ncor, nincor, nleft, nright)
                nival = int(nival * nival_fact + 0.5)

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
            if flags >= 2048:
                nsecondary += 1
                continue
            next_reference_start = int(toks[7])

            if '-' in coord_str:
                # paired end read, coordinate gives expected positions for each end
                c1, c2 = coord_str.split("-")

                if (flags & 3) != 3:
                    nnot_proper_pair += 1
                    continue  # not paired or not proper pair

                pos1, pos2 = int(c1), int(c2)

                # only use left end of reads, but check that right end is in
                # correct location
                if pos < next_reference_start:
                    correct_map = (pos1 == pos and pos2 == next_reference_start)
                    #print('%s: pos1 (%d) == pos (%d) + 1 and pos2 (%d) == next_reference_start (%d) + 1' % ('Correct' if correct_map else 'INCORRECT', pos1, pos, pos2, next_reference_start))
                    nleft += 1
                else:
                    nright += 1
                    continue  # this is right end of read
            else:
                correct_map = int(coord_str) == pos

            hist[orig_name][0 if correct_map else 1] += 1
            if correct_map:
                ncor += 1
            else:
                nincor += 1

    if not keep:
        os.remove(remapped_sam_fn)

    return hist, nsecondary, nnot_proper_pair, ncor, nincor, nleft, nright


# Output a file with one line per input alignment:
# - Original MAPQ
# - Predicted MAPQ
# - # derived reads that aligned correctly
# - # derived reads that aligned incorrectly
def tabulate(bam_fn, out_fn, hist):
    mapq_re = re.compile('Zm:[iZ]:([0-9]+)')
    tab = defaultdict(int)
    proc = subprocess.Popen(['samtools', 'view', '-h', bam_fn], stdout=subprocess.PIPE)
    for ln in proc.stdout:
        if ln[0] == '@':
            continue
        toks = ln.split('\t')
        qname = toks[0]
        flags = int(toks[1])
        if flags >= 2048:
            continue
        if qname in hist:
            mapq = int(toks[4])
            cor, incor = hist[qname]
            orig_mapq = mapq_re.search(ln)
            if orig_mapq is None:
                raise RuntimeError('Could not parse original mapq from this line: ' + ln)
            orig_mapq = int(orig_mapq.group(1))
            tab[(mapq, orig_mapq, cor, incor)] += 1
    ret = proc.wait()
    if ret != 0:
        raise RuntimeError('samtools returned %d' % ret)
    with open(out_fn, 'wb') as ofh:
        for k, v in tab.items():
            mapq, orig_mapq, cor, incor = k
            print(','.join(map(str, [mapq, orig_mapq, cor, incor, v])), file=ofh)


def add_args(parser):
    parser.add_argument('--bam', metavar='path', type=str, required=True,
                        help='Input BAM file; aligner and arguments are inferred from this')
    parser.add_argument('--fastq', metavar='path', type=str, required=True,
                        help='WASP "remap" FASTQ file')
    parser.add_argument('--fastq2', metavar='path', type=str, required=False,
                        help='WASP "remap" FASTQ file (end 2 for paired-end)')
    parser.add_argument('--threads', metavar='N', type=int, default=16,
                        help='Use N simultaneous threads when aligning remap FASTQs')
    parser.add_argument('--output', metavar='path', type=str, required=True,
                        help='Place output table here')
    parser.add_argument('--skip', action='store_const', const=True, default=False,
                        help='Skip over alignment')
    parser.add_argument('--keep', action='store_const', const=True, default=False,
                        help='Keep SAM file')


def go(args):
    if not args.skip:
        print('Getting arguments from BAM', file=sys.stderr)
        aligner, aligner_args = args_from_bam(args.bam)
        print('  Aligner arguments: ' + str(aligner_args), file=sys.stderr)
        print('Aligning overlap FASTQ', file=sys.stderr)
        ofn = args.output
        fastq1_fn = args.fastq
        if not os.path.exists(fastq1_fn):
            raise RuntimeError('No such FASTQ as "%s"' % fastq1_fn)
        fastq2_fn = args.fastq2
        if fastq2_fn is not None and not os.path.exists(fastq2_fn):
            raise RuntimeError('No such FASTQ as "%s"' % fastq2_fn)
        ofn += '.sam'
        ofn = join(os.path.dirname(ofn), '.' + os.path.basename(ofn))
        if aligner == 'SNAP':
            if not is_exe(snap_exe):
                raise RuntimeError('No snap exe at: "%s"' % snap_exe)
            remap_sam_fn = align_fastq_snap(fastq1_fn, fastq2_fn, aligner_args, args.threads, ofn)
        elif aligner == 'bowtie2':
            if not is_exe(bowtie2_exe):
                raise RuntimeError('No bowtie2 exe at: "%s"' % bowtie2_exe)
            remap_sam_fn = align_fastq_bowtie2(fastq1_fn, fastq2_fn, aligner_args, args.threads, ofn)
        elif aligner == 'bwa':
            remap_sam_fn = align_fastq_bwa(fastq1_fn, fastq2_fn, aligner_args, args.threads, ofn)
        else:
            if not is_exe(bwa_exe):
                raise RuntimeError('No bwa exe at: "%s"' % bwa_exe)
            raise RuntimeError('Unknown aligner ID: "%s"' % aligner)
        print('Scanning alignments', file=sys.stderr)
    else:
        remap_sam_fn = args.bam
    hist, nsecondary, nnot_proper_pair, ncor, nincor, nleft, nright = scan_remapped_bam(remap_sam_fn, keep=args.keep or args.skip)
    print_update(nsecondary, nnot_proper_pair, ncor, nincor, nleft, nright)
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
                    'correctly and incorrectly')

    add_args(_parser)
    _parser.add_argument('--profile', action='store_const', const=True, default=False, help='Print profiling info')
    _args = _parser.parse_args(sys.argv[1:])

    go_profile(_args)
