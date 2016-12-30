"""
You give me some paired-end reads and some SAM files from various
aligners having been run on those reads.  I match up the results of the
various aligners and output a new set of reads that are essentially a
subset of the input reads, modified in two key ways:


Iterate through all the SAM files in tandem.  Take one of the results round-robin and turn it into a string that Qtip will understand.

"""

from __future__ import print_function
import os
import sys
import gzip
from operator import itemgetter


def left_soft_clipped(cigar):
    """ Given a cigar string, return the indicated amount of soft
        clipping on the left side """
    n = 0
    for c in cigar:
        if str.isdigit(c):
            n *= 10
            n += (ord(c) - ord('0'))
        else:
            if c == 'S':
                return n
            return 0  # leftmost op was something other than S
    return 0


round_robin = 0


def handle_lines(sam_linetoks, fastq_lines, ofhs):
    """


    Example: 11_25006153_25006410_0:0:0_0:0:0_100_100_1_1/1
             ^ refid    ^ frag end (1-based) ^ len1  ^ flip
                ^ frag start (1-based)            ^ len2
                                  (bunch of 0s)

    fprintf(fpo[j], "@%s_%u_%u_%d:%d:%d_%d:%d:%d_%llx/%d\n", ks->name.s, ext_coor[0]+1, ext_coor[1]+1,
            n_err[0], n_sub[0], n_indel[0], n_err[1], n_sub[1], n_indel[1],
            (long long)ii, j==0? is_flip+1 : 2-is_flip);

    1. ks->name.s: ref name
    2. ext_coor[0]+1: leftmost
    3. ext_coor[1]+1: rightmost
    4. n_err[0], 5. n_sub[0], 6. n_indel[0]
    7. n_err[1], 8. n_sub[1], 9. n_indel[1]
    10. ii is just the # of the pair
    11. mate #
    """
    global round_robin

    sam_linetoks1 = sam_linetoks[0][round_robin]  # one mate
    sam_linetoks2 = sam_linetoks[1][round_robin]  # other mate from same aligner
    round_robin = (round_robin + 1) % len(sam_linetoks[0])  # pick next round robin

    flags1, flags2 = int(sam_linetoks1[1]), int(sam_linetoks2[1])
    assert (flags1 & 1) != 0 and (flags2 & 1) != 0
    if (flags2 & 64) != 0:
        sam_linetoks1, sam_linetoks2 = sam_linetoks2, sam_linetoks1

    refid = sam_linetoks1[2]
    assert sam_linetoks2[2] == refid

    pos1, pos2 = int(sam_linetoks1[3]), int(sam_linetoks2[3])
    len1, len2 = len(sam_linetoks1[10]), len(sam_linetoks2[10])
    cigar1, cigar2 = sam_linetoks1[5], sam_linetoks2[5]
    clip1, clip2 = left_soft_clipped(cigar1), left_soft_clipped(cigar2)
    pos1 -= clip1
    pos2 -= clip2
    flip = '0' if pos1 <= pos2 else '1'
    lpos = min(pos1, pos2)
    rpos = max(pos1 + len1 - 1, pos2 + len2 - 1)

    name = '@%s_%d_%d_0:0:0_0:0:0_%d_%d_%s_1' % (refid, lpos, rpos, len1, len2, flip)
    ofhs[0].write('%s/1\n%s+\n%s' % (name, fastq_lines[0][0], fastq_lines[0][1]))
    ofhs[1].write('%s/2\n%s+\n%s' % (name, fastq_lines[1][0], fastq_lines[1][1]))


def all_aligned(linetoks):
    """ Given list of lists, where outer list has alignments from
        different tools and inner has the line tokens for that tool,
        return true iff flags indicate all alignments are paired and
        concordant. """
    flags = map(int, map(itemgetter(1), linetoks))
    return all(map(lambda x: (x & 3) == 3 and (x & 4) == 0, flags))


def next_pair_fastq(fastq_fhs):
    """ Given list of filehandles for (hopefully parallel) FASTQ files,
        return list of pairs, where each pair is seq and qual from one
        file """
    fastq_lines = []
    first_nm = None
    for fastq_fh in fastq_fhs:
        nm = fastq_fh.readline()
        seq = fastq_fh.readline()
        nm2 = fastq_fh.readline()
        qual = fastq_fh.readline()
        if len(qual) == 0:
            raise RuntimeError('FASTQ ended before we expected')
        assert nm[0] == '@'
        assert nm2[0] == '+'
        if first_nm is None:
            first_nm = nm
        elif nm != first_nm and nm[:-2] != first_nm[:-2]:
            raise RuntimeError('Parallel FASTQ read names did not match:\n%s\n%s\n' % (nm, first_nm))
        fastq_lines.append((seq, qual))
    return fastq_lines


def next_pair_sam(sam_fhs):
    """ Return tokenized version of next non-secondary SAM record from several
        sam files """
    sam_lines = []
    for fh in sam_fhs:
        ln = '@'
        toks = None
        # skip headers and secondary alignments
        while ln[0] == '@' or int(toks[1]) >= 2048:
            ln = fh.readline()
            if len(ln) == 0:
                break
            toks = ln.split('\t')
        sam_lines.append(ln if len(ln) > 0 else None)
    if all(map(lambda x: x is None, sam_lines)):
        for fh in sam_fhs:
            fh.close()
        return True, None  # done
    if any(map(lambda x: x is None, sam_lines)):
        raise RuntimeError('Files finished out of sync')
    linetoks = list(map(lambda x: x.rstrip().split('\t'), sam_lines))
    for i in range(1, len(linetoks)):
        if linetoks[i][0] != linetoks[0][0]:
            raise RuntimeError('Read name mismatch between two SAM files: %s vs %s' %
                               (linetoks[i][0], linetoks[0][0]))
    return False, linetoks


def openex(fn, mode="rb"):
    if fn.endswith(".gz"):
        import gzip
        return gzip.open(fn, mode)
    else:
        return open(fn, mode)


def handle_sams(sam_fns, fastq_fns, output_fns):
    """
    Iterate through sams and FASTQs in lock step.
    """
    sam_fhs = map(openex, sam_fns)
    fastq_fhs = map(openex, fastq_fns)
    ln = None
    nreads = 0
    nreads_ival = 2
    nreads_ival_frac = 1.2
    with open(output_fns[0], 'wb') as ofh1:
        with open(output_fns[1], 'wb') as ofh2:
            while True:
                def close_all():
                    for fastq_fh in fastq_fhs:
                        fastq_fh.close()
                done, sam_linetoks1 = next_pair_sam(sam_fhs)
                if done:
                    close_all()
                    break
                done, sam_linetoks2 = next_pair_sam(sam_fhs)
                if done:
                    close_all()
                    break
                assert sam_linetoks1 is not None
                assert sam_linetoks2 is not None
                if ln is None:
                    ln = len(sam_linetoks1)
                fastq_lines = next_pair_fastq(fastq_fhs)
                assert fastq_lines is not None
                assert ln == len(sam_linetoks1)
                assert ln == len(sam_linetoks2)
                if all_aligned(sam_linetoks1) and all_aligned(sam_linetoks2):
                    # Write corresponding FASTQ records to output files
                    handle_lines([sam_linetoks1, sam_linetoks2],
                                 fastq_lines, [ofh1, ofh2])
                nreads += 1
                if nreads == nreads_ival:
                    print('Processed %d reads' % nreads)
                    nreads_ival = int(nreads_ival * nreads_ival_frac) + 1


def add_args(parser):
    parser.add_argument('--sam', metavar='path', type=str, nargs='+', required=True,
                        help='Input SAM file; aligner and arguments are inferred from this')
    parser.add_argument('--fastq1', metavar='path', type=str, required=True,
                        help='Input fastq1')
    parser.add_argument('--fastq2', metavar='path', type=str, required=True,
                        help='Input fastq2')
    parser.add_argument('--output1', metavar='path', type=str, required=True,
                        help='Place end-1 output FASTQ here')
    parser.add_argument('--output2', metavar='path', type=str, required=True,
                        help='Place end-2 output FASTQ here')


def go(args):
    for bm in args.sam:
        if not os.path.exists(bm):
            raise RuntimeError('Input SAM file "%s" does not exist' % bm)
    if not os.path.exists(args.fastq1):
        raise RuntimeError('Input FASTQ file "%s" does not exist' % args.fastq1)
    if not os.path.exists(args.fastq2):
        raise RuntimeError('Input FASTQ file "%s" does not exist' % args.fastq2)
    handle_sams(args.sam, [args.fastq1, args.fastq2], [args.output1, args.output2])


def go_profile(args):
    if args.profile:
        import cProfile
        cProfile.run('go(args)')
    else:
        go(args)


if __name__ == "__main__":

    import argparse

    _parser = argparse.ArgumentParser(
        description="""
Given SAM files from many different aligners, but all derived from
same input FASTQ, generate new pair of unpaired FASTQ files with only
reads aligned by all tools and with correctness info installed in
read name.  Later, Qtip can use that correctness info""")

    add_args(_parser)
    _parser.add_argument('--profile', action='store_const', const=True, default=False, help='Print profiling info')
    _args = _parser.parse_args(sys.argv[1:])

    go_profile(_args)
