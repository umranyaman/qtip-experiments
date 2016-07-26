"""
Check whether there's good concordance between training and test data
in various ways.
"""

from __future__ import print_function
from collections import defaultdict
import sys
import re
import logging


"""
First pass:
- Collect MAPQ distributions for tandem & input
- Collect ZT:Z distributions for tandem & input

Second pass:
- For tandem and input:
  + For each read, check for incorrectness.  If incorrect:
    - Record vector of percentiles for MAPQ and ZT:Zs

Finally:
- For tandem and input:
  + Sort the percentile vector by MAPQ percentile and record

Example usage:
pypy eval_concordance.py $HOME/final.sam $HOME/tandem_unp.sam inp.incor tan.incor

TODO:
- Are we doing paired-end correctly for the various name formats?
"""


def pass1_line(ln, mapq_dist, ztz_dist):
    """ Extracts MAPQ and ZT:Z fields; sticks them in histograms """
    if ln[0] == '@':
        return
    toks = ln.rstrip().split('\t')
    if toks[1] == '4':
        return
    mapq_dist[int(toks[4])] += 1  # MAPQ
    for tok in toks[12:]:
        if tok.startswith('ZT:Z:'):
            for i, ztok in enumerate(tok[6:].split(',')):
                ztz_dist[i][float(ztok)] += 1  # ZT:Z
            break


def pass1_fh(ifh, tfh):
    """ Histograms all the input and tandem alignments """
    mapq_dist_inp = defaultdict(int)
    ztz_dist_inp = defaultdict(lambda: defaultdict(int))
    for ln in ifh:
        pass1_line(ln, mapq_dist_inp, ztz_dist_inp)
    mapq_dist_tan = defaultdict(int)
    ztz_dist_tan = defaultdict(lambda: defaultdict(int))
    for ln in tfh:
        pass1_line(ln, mapq_dist_tan, ztz_dist_tan)
    return mapq_dist_inp, ztz_dist_inp, mapq_dist_tan, ztz_dist_tan


def pass1_fn(input_fn, tandem_fn):
    """ Histograms all the input and tandem alignments """
    with open(input_fn) as ifh:
        with open(tandem_fn) as tfh:
            return pass1_fh(ifh, tfh)

"""
Example: 10_26049747_26049846_0:0:0_0:0:0_100_100_0_3999999
offset is 0-based?
"""
_wgsimex_re = re.compile('(.+)_([^_]+)_([^_]+)_([^:]+):([^:]+):([^_]+)_([^:]+):([^:]+):([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^/]+).*')
#                           1   2       3       4       5       6       7       8       9       10      11      12      13


def name_is_extended_wgsim(nm):
    return _wgsimex_re.match(nm) is not None


def pos_from_extended_wgsim(name, mate2=False):
    res = _wgsimex_re.match(name)
    refid, fragst1, fragen1 = res.group(1), int(res.group(2))-1, int(res.group(3))-1
    len1, len2 = int(res.group(10)), int(res.group(11))
    flip = res.group(12) == '1'
    ln = len2 if mate2 else len1
    if flip == mate2:
        return refid, fragst1, True
    else:
        return refid, fragen1 - (ln-1), False

"""
Example: qsim!:GL000229.1:+:6005:100:u
pretty sure offset is 0-based
"""
_qsim_re = re.compile('qsim!:([^:]+):([+-]):([^:]+):.*')


def name_is_qsim(nm):
    return _qsim_re.match(nm) is not None


def pos_from_qsim(name):
    res = _qsim_re.match(name)
    return res.group(1), int(res.group(3)), res.group(2) == '+'


"""
Example: !h!chr9!118085975!+!50!0
offset is 0-based
"""
_hint_re = re.compile('!h!([^!]+)!([0-9]+)!([+-])!([0-9]+)!([0-9]+)')


def name_is_hint(name):
    return _hint_re.match(name) is not None


def pos_from_hint(name):
    res = _hint_re.match(name)
    return res.group(1), int(res.group(2)), res.group(3) == '+'


"""
Example:
hg38_50nt_mason1_unp.fastq.000999999 contig=chr9 haplotype=1 length=50 orig_begin=118085976 orig_end=118086026 snps=0 indels=0 haplotype_infix=ATGACTCTTGAAGCTGGGCGCAGTGGCTCATGCCTGTAATCCTAGCACTT edit_string=MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM strand=forward
"""
_mason_re = re.compile('[^ ]+ contig=([^ ]+) .* orig_begin=([^ ]+) .* strand=([fr]).*')


def name_is_mason1(name):
    return _mason_re.match(name) is not None


def pos_from_mason1(name):
    res = _hint_re.match(name)
    return res.group(1), int(res.group(2)), res.group(3) == 'f'


def same_pos(pos1, pos2, wiggle=30):
    """ Returns true when the two positions are basically the same """
    refid1, pos1, strand1 = pos1
    refid2, pos2, strand2 = pos2
    if refid1 != refid2 or strand1 != strand2:
        return False
    return abs(pos1 - pos2) < wiggle


def is_correct(toks, wiggle=30):
    """ Checks whether alignment, tokenized in toks, is correct """
    flags = int(toks[1])
    aligned_pos = (toks[2], int(toks[3])-1, (flags & 16) == 0)
    paired = (flags & 1) != 0
    mate2 = paired and (flags & 128) != 0
    if name_is_extended_wgsim(toks[0]):
        true_pos = pos_from_extended_wgsim(toks[0], mate2)
    elif name_is_qsim(toks[0]):
        true_pos = pos_from_qsim(toks[0])
    elif name_is_mason1(toks[0]):
        true_pos = pos_from_mason1(toks[0])
    elif name_is_hint(toks[0]):
        true_pos = pos_from_hint(toks[0])
    else:
        raise RuntimeError('Name was not formatted as expected: "%s"' % toks[0])
    return same_pos(true_pos, aligned_pos, wiggle=wiggle)


def percentile_in(item, dist):
    # get item's percentile in empirical distribution dist
    pass


def percentileize(dist):
    last_v = 0
    cum = {}
    for k, v in sorted(dist.items()):
        assert k not in cum
        cum[k] = last_v
        last_v += v
    pct_dict = {}
    for k in dist.keys():
        pct_dict[k] = (cum[k], float(cum[k])/last_v, dist[k])
    return pct_dict


def pass2_fh(fh, mapq_dist, ztz_dist, ofh):
    logging.info('  Calculating percentiles')
    mapq_pctile_dict = percentileize(mapq_dist)
    ztz_pctile_dict = {x: percentileize(v) for x, v in ztz_dist.items()}
    incors = []
    for ln in fh:
        if ln[0] == '@':
            continue
        toks = ln.rstrip().split('\t')
        if toks[1] == '4':
            continue
        if not is_correct(toks):
            mapq = int(toks[4])
            assert mapq in mapq_pctile_dict
            mapq_pctile = mapq_pctile_dict[mapq]
            ls = [mapq_pctile[1]]
            for tok in toks[12:]:
                if tok.startswith('ZT:Z:'):
                    for i, ztok in enumerate(tok[6:].split(',')):
                        flztok = float(ztok)
                        assert flztok in ztz_pctile_dict[i]
                        ztz_pctile = ztz_pctile_dict[i][flztok]
                        ls.append(ztz_pctile[1])
                    break
            incors.append(ls)
    logging.info('  Sorting and printing')
    for ls in sorted(incors, reverse=True):
        ofh.write(','.join(map(str, ls)) + '\n')


def pass2_fn(input_fn, tandem_fn, mapq_dist_inp, ztz_dist_inp, mapq_dist_tan, ztz_dist_tan, ofn_inp, ofn_tan):
    with open(ofn_inp, 'w') as ofh:
        with open(input_fn) as ifh:
            pass2_fh(ifh, mapq_dist_inp, ztz_dist_inp, ofh)
    with open(ofn_tan, 'w') as ofh:
        with open(tandem_fn) as tfh:
            pass2_fh(tfh, mapq_dist_tan, ztz_dist_tan, ofh)


def go():
    format_str = '%(asctime)s:%(levelname)s:%(message)s'
    level = logging.INFO
    logging.basicConfig(format=format_str, datefmt='%m/%d/%y-%H:%M:%S', level=level)

    logging.info('Pass 1')
    mapq_dist_inp, ztz_dist_inp, mapq_dist_tan, ztz_dist_tan = pass1_fn(sys.argv[1], sys.argv[2])
    logging.info('Pass 2')
    pass2_fn(sys.argv[1], sys.argv[2], mapq_dist_inp, ztz_dist_inp, mapq_dist_tan, ztz_dist_tan, sys.argv[3], sys.argv[4])
    logging.info('Done')


if __name__ == '__main__':
    go()
