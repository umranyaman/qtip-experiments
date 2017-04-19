"""
cat sam_file.sam | evaluate.py {human|mouse}

Iterate through a SAM file produced by the Makefile and count:

- # category-1a errors (contaminant aligned to primary reference)
- # category-1b errors (same-species non-reference sequence aligned to primary reference)
- # category-2 errors (reference-derived reads that failed to align)
- # category-3 errors (aligned to reference but not to point of origin)
- correct alignments of both kinds (target and contamination)
"""

from __future__ import print_function
import sys
import re


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


if len(sys.argv) <= 1:
    raise RuntimeError('Specify either "human" or "mouse" as argument')


human_chrs = {'1', '2', '3', '4', '5', '6', '7', '8', '9',
              '10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
              '20', '21', '22', 'MT', 'X', 'Y'}


mouse_chrs = {'1', '2', '3', '4', '5', '6', '7', '8', '9',
              '10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
              'MT', 'X', 'Y'}


def is_human_chr(st):
    return st in human_chrs or st[:2] == 'KI' or st[:2] == 'GL'


def is_mouse_chr(st):
    return st in mouse_chrs or st[:2] == 'JH' or st[:2] == 'GL'


wiggle = 30
cat1a, cat1b, cat2, cat3 = 0, 0, 0, 0
cor_target = 0  # # target reads correctly aligned to target
cor_contam = 0
for ln in sys.stdin:
    if ln[0] == '@':
        continue
    toks = ln.split('\t')
    from_contaminant = toks[0][0] != 'r'
    from_chm = toks[0].startswith('utg718000') or toks[0].startswith('JSAF020')
    assert not from_contaminant or not from_chm
    unal = (int(toks[1]) & 4) != 0
    if unal and not from_contaminant and not from_chm:
        cat2 += 1  # incorrectly failed to align to target
    elif unal:
        cor_contam += 1  # correctly failed to align to target
    else:
        assert sys.argv[1] != 'human' or is_human_chr(toks[2])
        assert sys.argv[1] != 'mouse' or is_mouse_chr(toks[2])
        if from_contaminant:
            cat1a += 1  # incorrectly aligned to target
        elif from_chm:
            cat1b += 1  # incorrectly aligned to target
        elif not is_correct(toks, wiggle):
            cat3 += 1  # correctly aligned to target, but to wrong locus
        else:
            cor_target += 1  # correct

cat1 = cat1a + cat1b
err = cat1 + cat2 + cat3
tot = err + cor_target + cor_contam

print('total=%d' % tot, file=sys.stderr)
print('error=%d, %0.04f%%' % (err, float(100*err)/tot), file=sys.stderr)

print('type,count,pct_total,pct_error')
print('1,%d,%0.04f,%0.04f' % (cat1, float(100*cat1)/tot, float(100*cat1)/err))
print('1a,%d,%0.04f,%0.04f' % (cat1a, float(100*cat1a)/tot, float(100*cat1a)/err))
print('1b,%d,%0.04f,%0.04f' % (cat1b, float(100*cat1b)/tot, float(100*cat1b)/err))
print('2,%d,%0.04f,%0.04f' % (cat2, float(100*cat2)/tot, float(100*cat2)/err))
print('3,%d,%0.04f,%0.04f' % (cat3, float(100*cat3)/tot, float(100*cat3)/err))
