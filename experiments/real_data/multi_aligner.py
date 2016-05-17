from __future__ import print_function

import sys
import re
import os
import numpy
import time
from collections import defaultdict

__author__ = 'langmead'

"""
Ultimately, this tool helps us compare MAPQ assignments for real-word reads
and for many tools by asking how well the MAPQs reflect a proxy for the
"truth."  The proxy is the degree to which the various tools agree on the
reported alignment.

We consider the aligners/SAM files to be in "tiers."  E.g. maybe the SAMs
correspond to various levels of sensitivity for the same tool, so we care more
about whether an alignment is in agreement with the "higher sensitivity" tiers
than with the other tiers.

Given SAM files from multiple alignment tools, sort them by read id, iterate
through all in tandem and, for each read:

(a) if any of the tools failed to alignment the read, in which case, skip
(b) form a matrix of which tools agree with each other
(c) for each row (aligner), ask how many times the aligner agrees with a
    higher-tier aligner and how many times it agrees with a same-tier aligner

For every alignment we have information like this

Tool   Tier   MAPQ  Row of similarity matrix    Agreement w/  Agreement w/  Agreement w/
                                                higher tiers    own tier       others
bt2_vs    1     11   True,  True, False, False       0             1             2
bt2_s     2     10   True,  True, False, False       1             1             2
bt2_f     3     30  False, False,  True, False       0             1             1
bt2_v     4     30  False, False, False,  True       0             1             1

We then output several ROC tables = 6 * # tools.  6 because we produce 1 for
both "strict" and "loose" definitions of correctness and using (a) old, (b)
new, and (c) new, rounded MAPQ values.

"""

VERSION = '0.0.0'


def add_args(parser):
    # Overall arguments
    parser.add_argument('--wiggle', metavar='int', type=int, required=False, default=10,
                        help='# bases difference there can be between two alignments that are the same')
    parser.add_argument('--sam', metavar='paths', type=str, nargs='+', required=True,
                        help='SAM file(s) containing aligned reads')
    parser.add_argument('--name', metavar='strings', type=str, nargs='+', required=True,
                        help='Tool names')
    parser.add_argument('--prefix', type=str, required=False, help='Prefix for output files')
    parser.add_argument('--suffix', type=str, required=False, help='Suffix for output files')
    parser.add_argument('--tier', metavar='ints', type=int, nargs='+', required=True,
                        help='Tiers, lower being more accurate (in expectation)')

mapq_re = re.compile('Zm:[iZ]:([0-9]+)')
mapq_prec_re = re.compile('Zp:Z:([.0-9]+)')
digre = re.compile('([01-9]+)')


def cigar_to_list(cigar_string):
    spl = digre.split(cigar_string)
    return zip(spl[2::2], map(int, spl[1::2]))


def parse_sam_loc_mapq(st):
    toks = st.split('\t')
    flags = int(toks[1])
    if (flags & 4) != 0:
        return None  # failed to align
    rname = toks[2]
    rpos = int(toks[3])
    clist = cigar_to_list(toks[5])
    lclip = 0 if clist[0][0] != 'S' else clist[0][1]
    mapq = int(toks[4])  # predicted mapq
    # now need to go get the original mapq
    mapq_re_ma = mapq_re.search(st)
    if mapq_re_ma is None:
        raise RuntimeError('')
    mapq_orig = mapq_re_ma.group(1)
    mapq_pred_re_ma = mapq_prec_re.search(st)
    if mapq_pred_re_ma is None:
        raise RuntimeError('')
    mapq_prec = mapq_pred_re_ma.group(1)
    return rname, rpos - lclip, int(mapq), int(mapq_orig), float(mapq_prec)


def compare_2_alns(a, b, wiggle=50):
    """ Return true iff the two alignments seem to be to the same place """
    arefid, apos, _, _, _ = a
    brefid, bpos, _, _, _ = b
    return arefid == brefid and abs(apos - bpos) < wiggle


def same_alns_count(ls, wiggle=50):
    """ Return a list of ints indicating the number of tools "agreeing" with
        this alignment. """
    sim = [1] * len(ls)
    for i in range(len(ls)):
        for j in range(i+1, len(ls)):
            if compare_2_alns(ls[i], ls[j], wiggle=wiggle):
                sim[i] += 1
                sim[j] += 1
    return sim


def same_alns_matrix(ls, wiggle=50):
    """ Return a numpy matrix of bools indicating which tools agree with each
        other. """
    sim = numpy.zeros((len(ls), len(ls)), dtype=bool)
    for i in range(len(ls)):
        sim[i, i] = True
        for j in range(i+1, len(ls)):
            sim[i, j] = sim[j, i] = compare_2_alns(ls[i], ls[j], wiggle=wiggle)
    return sim


def same_alns_tiered(ls, better_tiers, equal_tiers, wiggle=50):
    """ Return a numpy matrix of bools indicating which tools agree with each
        other. """
    sim = numpy.zeros((len(ls), len(ls)), dtype=bool)
    tier_results = []

    for i in range(len(ls)):
        sim[i, i] = True
        for j in range(i+1, len(ls)):
            sim[i, j] = sim[j, i] = compare_2_alns(ls[i], ls[j], wiggle=wiggle)

    for i in range(len(ls)):
        num_better = numpy.sum(sim[i, better_tiers[i]])
        num_equal = numpy.sum(sim[i, equal_tiers[i]])
        num_oall = numpy.sum(sim[i])
        correct_l = num_better == len(better_tiers[i])
        correct_s = correct_l and num_equal == len(equal_tiers[i])
        tier_results.append((num_better, num_equal, num_oall, correct_l, correct_s))

    return tier_results


def preprocess_tiers(tiers):
    better_tiers, equal_tiers = [], []
    for i in range(len(tiers)):
        better_tier = []
        equal_tier = []
        for j in range(len(tiers)):
            if tiers[j] < tiers[i]:
                better_tier.append(j)
            elif tiers[j] == tiers[i]:
                equal_tier.append(j)
        better_tiers.append(better_tier)
        equal_tiers.append(equal_tier)
    return better_tiers, equal_tiers


def run(cmd):
    print(cmd, file=sys.stderr)
    return os.system(cmd)


def preprocess_sm(fn):
    print("Preprocessing %s..." % fn, file=sys.stderr)
    run("which samtools")
    if not os.path.exists(fn + '.bam'):
        ret = run("samtools view -b %s > %s" % (fn, fn + '.bam'))
        if ret != 0:
            raise RuntimeError('samtools view-to-bam failed')
        os.remove(fn)
    if not os.path.exists(fn + '.sorted.bam'):
        ret = run("samtools sort -n %s %s" % (fn + '.bam', fn + '.sorted'))
        if ret != 0:
            raise RuntimeError('samtools sort failed')
        os.remove(fn + '.bam')
    if not os.path.exists(fn + '.sorted.sam'):
        ret = run("samtools view %s > %s" % (fn + '.sorted.bam', fn + '.sorted.sam'))
        if ret != 0:
            raise RuntimeError('samtools view-to-sam failed')
        os.remove(fn + '.sorted.bam')
    print("DONE", file=sys.stderr)
    return fn + '.sorted.sam'


def write_tally_to_roc(tally, fn):
    """ Return the ROC table given a list of pcors and a parallel list of
        correct/incorrect booleans. """
    cum_cor, cum_incor = 0, 0
    with open(fn, 'wb') as fh:
        fh.write(','.join(['mapq', 'cor', 'incor', 'cum_cor', 'cum_incor']) + '\n')
        for p in sorted(tally.keys(), reverse=True):
            ncor, nincor = tally[p]
            cum_cor += ncor
            cum_incor += nincor
            fh.write(','.join(map(str, [p, ncor, nincor, cum_cor, cum_incor])) + '\n')


def auc(tally):
    """ Calculate area under curve given predictions and correct/incorrect
        information. """
    area, tot_cor, tot_incor = 0, 0, 0
    last_tot_cor, last_tot_incor = 0, 0
    for pcor, ci in sorted(tally.items(), reverse=True):
        tot_cor += ci[0]
        tot_incor += ci[1]
        cor_diff = tot_cor - last_tot_cor
        incor_diff = tot_incor - last_tot_incor
        if incor_diff > 0:
            area += (0.5 * cor_diff * incor_diff)
            area += last_tot_cor * incor_diff
        last_tot_cor, last_tot_incor = tot_cor, tot_incor
    return area


def ranking_error(tally):
    """ Return the ranking error given a list of pcors and a parallel list of
        correct/incorrect booleans.  Round off to nearest MAPQ first if
        rounded=True.  """
    err, sofar = 0, 0
    # from low-confidence to high confidence
    for p, tup in sorted(tally.items()):
        ncor, nincor = tup
        ntot = ncor + nincor
        assert ntot > 0
        if nincor > 0:
            # spread error over this grouping of tied pcors
            frac = float(nincor) / ntot
            assert frac <= 1.0
            err += frac * sum(range(sofar, sofar + ntot))
        sofar += ntot
    return err


def go(args):

    def decorate_filename(fn):
        return ('' if args['prefix'] is None else args['prefix']) + fn + ('' if args['suffix'] is None else args['suffix'])

    names = args['name']
    sams = map(lambda x: open(x, 'rb'), map(preprocess_sm, args['sam']))
    better_tiers, equal_tiers = preprocess_tiers(map(int, args['tier']))
    assert len(sams) == len(better_tiers)
    assert len(sams) == len(equal_tiers)
    assert len(sams) == len(names)
    rocs = {}
    for nm in args['name']:
        rocs[nm + "_loose_orig"] = defaultdict(lambda: [0, 0])
        rocs[nm + "_loose_int"] = defaultdict(lambda: [0, 0])
        rocs[nm + "_loose_dec"] = defaultdict(lambda: [0, 0])
        rocs[nm + "_strict_orig"] = defaultdict(lambda: [0, 0])
        rocs[nm + "_strict_int"] = defaultdict(lambda: [0, 0])
        rocs[nm + "_strict_dec"] = defaultdict(lambda: [0, 0])
    ival = 200000
    lni = 0
    itime = time.time()
    while True:
        lns, rdnames = [], []
        for sam in sams:
            while True:
                ln = sam.readline()
                if len(ln) == 0:
                    lns.append(None)
                    break
                elif ln[0] == '@':
                    continue
                toks = ln.split('\t')
                flags = int(toks[1])
                if (flags & 256) != 0:  # skip secondary
                    continue
                if (flags & 2048) != 0:  # skip spliced
                    continue
                lns.append(ln)
                rdnames.append(toks[0])
                break

        if not all(map(lambda x: rdnames[0] == x, rdnames)):
            print("Warning: read names don't match: %s" % str(rdnames), file=sys.stderr)
        assert len(lns) == len(sams)
        if all(map(lambda x: x is None, lns)):
            break  # reached end of all the files
        # if we reach end of one, should reach end of all
        assert not any(map(lambda x: x is None, lns))
        parsed_sam = map(parse_sam_loc_mapq, lns)
        if any(map(lambda x: x is None, parsed_sam)):
            continue  # at least one tool failed to align the read
        sim = same_alns_tiered(parsed_sam, better_tiers, equal_tiers, args['wiggle'])
        for i, tup in enumerate(zip(parsed_sam, names)):
            ps, nm = tup
            rname, pos, mapq, mapq_orig, mapq_prec = ps
            print (','.join(map(str, [i] + list(sim[i]) + [mapq, mapq_orig, mapq_prec])))
            correct_l, correct_s = sim[i][-2], sim[i][-1]
            rocs[nm + '_loose_orig'][mapq_orig][0 if correct_l else 1] += 1
            rocs[nm + '_strict_orig'][mapq_orig][0 if correct_s else 1] += 1
            rocs[nm + '_loose_int'][mapq][0 if correct_l else 1] += 1
            rocs[nm + '_strict_int'][mapq][0 if correct_s else 1] += 1
            rocs[nm + '_loose_dec'][mapq_prec][0 if correct_l else 1] += 1
            rocs[nm + '_strict_dec'][mapq_prec][0 if correct_s else 1] += 1
        lni += 1
        if args['upto'] is not None and lni >= args['upto']:
            break
        if lni % ival == 0:
            elapsed_time = time.time() - itime
            print('  Handled %d alignments across %d files (~%d per sec)...' %
                  (lni, len(sams), int(lni / elapsed_time)), file=sys.stderr)

    for k, l in rocs.items():
        fn = decorate_filename(k) + '.roc.csv'
        write_tally_to_roc(l, fn)
        print('Wrote ROC table output to "%s"' % fn, file=sys.stderr)

        fn = decorate_filename(k) + '.summ.csv'
        with open(fn, 'wb') as fh:
            fh.write('%f,%f\n' % (auc(l), ranking_error(l)))
        print('Wrote summary output to "%s"' % fn, file=sys.stderr)


def go_profile(args):
    if args['profile']:
        import cProfile
        cProfile.run('go(args)')
    else:
        go(args)

if __name__ == "__main__":

    import argparse

    _parser = argparse.ArgumentParser(
        description='Take a set of SAM files output by different read '
                    'alignment tools running on the same input reads.  The SAM'
                    'files should come from the Qsim system, and should '
                    'include the predicted and old MAPQs.')

    if '--version' in sys.argv:
        print('Tandem simulator, version ' + VERSION, file=sys.stderr)
        sys.exit(0)

    if '--test' in sys.argv:
        import unittest

        class MultiAlignerTests(unittest.TestCase):

            def test_parse_sam_loc_mapq_1(self):
                rname, pos, mapq, mapq_orig, mapq_prec = parse_sam_loc_mapq(
                    '''ERR050082.25000	161	1	1670466	4	100M	=	1608011	0	ATAACCGTGGTGGACAGCATCTGCACCGCACCTGCGGGAGGGAGGGGGCCGAAGACAAGAGGGAGAATCACCCCTCCCGTGCCTGCAGTGGGCGCCACAC	:EDFEHGFJA;LJGJKBKHNDKBJG?IIHHJC9ILLFNLIGG9BA9GKHGGK9C7IKHG@KHJ9=?D8IGI=G>?JFH><=;A;JICG<KG>D&1@9@7B	AS:i:-6	XS:i:-10	XN:i:0	XM:i:2	XO:i:	XG:i:0	NM:i:2	MD:Z:93T4C1	YT:Z:UP	Zm:Z:6	Zp:Z:4.399''')
                self.assertEqual(rname, '1')
                self.assertEqual(pos, 1670466)
                self.assertEqual(mapq, 4)
                self.assertEqual(mapq_orig, 6)
                self.assertAlmostEqual(mapq_prec, 4.399)

            def test_compare_2_alns_1(self):
                bt2_pysam = parse_sam_loc_mapq(
                    '''ERR050082.25000	161	1	1670466	4	100M	=	1608011	0	ATAACCGTGGTGGACAGCATCTGCACCGCACCTGCGGGAGGGAGGGGGCCGAAGACAAGAGGGAGAATCACCCCTCCCGTGCCTGCAGTGGGCGCCACAC	:EDFEHGFJA;LJGJKBKHNDKBJG?IIHHJC9ILLFNLIGG9BA9GKHGGK9C7IKHG@KHJ9=?D8IGI=G>?JFH><=;A;JICG<KG>D&1@9@7B	AS:i:-6	XS:i:-10	XN:i:0	XM:i:2	XO:i:	XG:i:0	NM:i:2	MD:Z:93T4C1	YT:Z:UP	Zm:Z:6	Zp:Z:4.399''')
                bwa_pysam = parse_sam_loc_mapq(
                    '''ERR050082.25000	163	1	1670466	41	100M	=	1670869	503	ATAACCGTGGTGGACAGCATCTGCACCGCACCTGCGGGAGGGAGGGGGCCGAAGACAAGAGGGAGAATCACCCCTCCCGTGCCTGCAGTGGGCGCCACAC	:EDFEHGFJA;LJGJKBKHNDKBJG?IIHHJC9ILLFNLIGG9BA9GKHGGK9C7IKHG@KHJ9=?D8IGI=G>?JFH><=;A;JICG<KG>D&1@9@7B	NM:i:2	MD:Z:93T4C1	AS:i:93	XS:i:88	XA:Z:1,+1607608,100M,3;	Zm:Z:27	Zp:Z:40.752''')
                self.assertTrue(compare_2_alns(bt2_pysam, bwa_pysam, wiggle=10))

            def test_compare_2_alns_2(self):
                bt2_pysam = parse_sam_loc_mapq(
                    '''ERR050082.24499	163	1	1639199	5	100M	=	1639570	471	CGGCTGAGACAGAGCCCGGATGCTGAGCTGGGAGGAGGCGTCGGGTGTCATGTGGGGGACAAGCCCACATCCACGTCCACCAGGCTGAGGAAATAACCTA	:FDGEHFFCAKIJGCKMJJHJDGIIKIJJHJILILLLNLIGJLKKGMHMJGK@J9IIJGG+LIIGJJJE?8FFIIGK+DCJIDKB4GDFF1*9A@C@*0+	AS:i:-5	XS:i:-5	XN:i:0	XM:i:2	XO:i:0	XG:i:	NM:i:2	MD:Z:91C7C0	YS:i:0	YN:i:-60	Yn:i:0	ZN:i:-60	Zn:i:0	YT:Z:CP	Zm:Z:4	Zp:Z:5.088''')
                bwa_pysam = parse_sam_loc_mapq(
                    '''ERR050082.24501	83	1	1639728	41	100M	=	1639307	-521	GGTGCTGCCACAGGCAGGATGCGGGCTCCGCTTCAGCTAAGCAACAAGTGTTCCCAAGAATGGATATGGAGGCTGGGCGCGGTGGCTCACGCCTGTAATC	&=7@@9CEC?I:KGFDFBJ?GBJHIK@JKEHKJHJLKI?ELIGHKIJIAIKJFLKGJJEJGKLKEGMNMNKLNLLJJILCNKMMIKKIDJIKKJKHFHBD	NM:i:1	MD:Z:0A99	AS:i:99	XS:i:89	XA:Z:1,-1576518,100M,3;	Zm:Z:27	Zp:Z:40.752''')
                self.assertTrue(not compare_2_alns(bt2_pysam, bwa_pysam, wiggle=10))
                self.assertTrue(compare_2_alns(bt2_pysam, bwa_pysam, wiggle=1000))

            def test_count_same_alns_1(self):
                a = parse_sam_loc_mapq(
                    '''ERR050082.24499	163	1	1639199	5	100M	=	1639570	471	CGGCTGAGACAGAGCCCGGATGCTGAGCTGGGAGGAGGCGTCGGGTGTCATGTGGGGGACAAGCCCACATCCACGTCCACCAGGCTGAGGAAATAACCTA	:FDGEHFFCAKIJGCKMJJHJDGIIKIJJHJILILLLNLIGJLKKGMHMJGK@J9IIJGG+LIIGJJJE?8FFIIGK+DCJIDKB4GDFF1*9A@C@*0+	AS:i:-5	XS:i:-5	XN:i:0	XM:i:2	XO:i:0	XG:i:	NM:i:2	MD:Z:91C7C0	YS:i:0	YN:i:-60	Yn:i:0	ZN:i:-60	Zn:i:0	YT:Z:CP	Zm:Z:4	Zp:Z:5.088''')
                b = parse_sam_loc_mapq(
                    '''ERR050082.24501	83	1	1639728	41	100M	=	1639307	-521	GGTGCTGCCACAGGCAGGATGCGGGCTCCGCTTCAGCTAAGCAACAAGTGTTCCCAAGAATGGATATGGAGGCTGGGCGCGGTGGCTCACGCCTGTAATC	&=7@@9CEC?I:KGFDFBJ?GBJHIK@JKEHKJHJLKI?ELIGHKIJIAIKJFLKGJJEJGKLKEGMNMNKLNLLJJILCNKMMIKKIDJIKKJKHFHBD	NM:i:1	MD:Z:0A99	AS:i:99	XS:i:89	XA:Z:1,-1576518,100M,3;	Zm:Z:27	Zp:Z:40.752''')
                c = parse_sam_loc_mapq(
                    '''ERR050082.24501	83	1	1639728	41	100M	=	1639307	-521	GGTGCTGCCACAGGCAGGATGCGGGCTCCGCTTCAGCTAAGCAACAAGTGTTCCCAAGAATGGATATGGAGGCTGGGCGCGGTGGCTCACGCCTGTAATC	&=7@@9CEC?I:KGFDFBJ?GBJHIK@JKEHKJHJLKI?ELIGHKIJIAIKJFLKGJJEJGKLKEGMNMNKLNLLJJILCNKMMIKKIDJIKKJKHFHBD	NM:i:1	MD:Z:0A99	AS:i:99	XS:i:89	XA:Z:1,-1576518,100M,3;	Zm:Z:27	Zp:Z:40.752''')
                sim = same_alns_count([a, b, c], wiggle=10)
                self.assertEqual([1, 2, 2], sim)
                sim = same_alns_count([a, b, c], wiggle=1000)
                self.assertEqual([3, 3, 3], sim)

            def test_count_same_alns_2(self):
                a = parse_sam_loc_mapq(
                    '''ERR050082.24499	163	1	1639199	5	100M	=	1639570	471	CGGCTGAGACAGAGCCCGGATGCTGAGCTGGGAGGAGGCGTCGGGTGTCATGTGGGGGACAAGCCCACATCCACGTCCACCAGGCTGAGGAAATAACCTA	:FDGEHFFCAKIJGCKMJJHJDGIIKIJJHJILILLLNLIGJLKKGMHMJGK@J9IIJGG+LIIGJJJE?8FFIIGK+DCJIDKB4GDFF1*9A@C@*0+	AS:i:-5	XS:i:-5	XN:i:0	XM:i:2	XO:i:0	XG:i:	NM:i:2	MD:Z:91C7C0	YS:i:0	YN:i:-60	Yn:i:0	ZN:i:-60	Zn:i:0	YT:Z:CP	Zm:Z:4	Zp:Z:5.088''')
                b = parse_sam_loc_mapq(
                    '''ERR050082.24501	83	1	1639728	41	100M	=	1639307	-521	GGTGCTGCCACAGGCAGGATGCGGGCTCCGCTTCAGCTAAGCAACAAGTGTTCCCAAGAATGGATATGGAGGCTGGGCGCGGTGGCTCACGCCTGTAATC	&=7@@9CEC?I:KGFDFBJ?GBJHIK@JKEHKJHJLKI?ELIGHKIJIAIKJFLKGJJEJGKLKEGMNMNKLNLLJJILCNKMMIKKIDJIKKJKHFHBD	NM:i:1	MD:Z:0A99	AS:i:99	XS:i:89	XA:Z:1,-1576518,100M,3;	Zm:Z:27	Zp:Z:40.752''')
                c = parse_sam_loc_mapq(
                    '''ERR050082.24501	83	1	1639728	41	100M	=	1639307	-521	GGTGCTGCCACAGGCAGGATGCGGGCTCCGCTTCAGCTAAGCAACAAGTGTTCCCAAGAATGGATATGGAGGCTGGGCGCGGTGGCTCACGCCTGTAATC	&=7@@9CEC?I:KGFDFBJ?GBJHIK@JKEHKJHJLKI?ELIGHKIJIAIKJFLKGJJEJGKLKEGMNMNKLNLLJJILCNKMMIKKIDJIKKJKHFHBD	NM:i:1	MD:Z:0A99	AS:i:99	XS:i:89	XA:Z:1,-1576518,100M,3;	Zm:Z:27	Zp:Z:40.752''')
                d = parse_sam_loc_mapq(
                    '''ERR050082.24501	83	1	1649728	41	100M	=	1639307	-521	GGTGCTGCCACAGGCAGGATGCGGGCTCCGCTTCAGCTAAGCAACAAGTGTTCCCAAGAATGGATATGGAGGCTGGGCGCGGTGGCTCACGCCTGTAATC	&=7@@9CEC?I:KGFDFBJ?GBJHIK@JKEHKJHJLKI?ELIGHKIJIAIKJFLKGJJEJGKLKEGMNMNKLNLLJJILCNKMMIKKIDJIKKJKHFHBD	NM:i:1	MD:Z:0A99	AS:i:99	XS:i:89	XA:Z:1,-1576518,100M,3;	Zm:Z:27	Zp:Z:40.752''')
                e = parse_sam_loc_mapq(
                    '''ERR050082.24501	83	1	1649728	41	100M	=	1639307	-521	GGTGCTGCCACAGGCAGGATGCGGGCTCCGCTTCAGCTAAGCAACAAGTGTTCCCAAGAATGGATATGGAGGCTGGGCGCGGTGGCTCACGCCTGTAATC	&=7@@9CEC?I:KGFDFBJ?GBJHIK@JKEHKJHJLKI?ELIGHKIJIAIKJFLKGJJEJGKLKEGMNMNKLNLLJJILCNKMMIKKIDJIKKJKHFHBD	NM:i:1	MD:Z:0A99	AS:i:99	XS:i:89	XA:Z:1,-1576518,100M,3;	Zm:Z:27	Zp:Z:40.752''')
                sim = same_alns_count([a, b, c, d, e], wiggle=10)
                self.assertEqual([1, 2, 2, 2, 2], sim)
                sim = same_alns_count([a, b, c, d, e], wiggle=1000)
                self.assertEqual([3, 3, 3, 2, 2], sim)
                sim = same_alns_count([a, b, c, d, e], wiggle=20000)
                self.assertEqual([5, 5, 5, 5, 5], sim)

            def test_same_alns_matrix_1(self):
                a = parse_sam_loc_mapq(
                    '''ERR050082.24499	163	1	1639199	5	100M	=	1639570	471	CGGCTGAGACAGAGCCCGGATGCTGAGCTGGGAGGAGGCGTCGGGTGTCATGTGGGGGACAAGCCCACATCCACGTCCACCAGGCTGAGGAAATAACCTA	:FDGEHFFCAKIJGCKMJJHJDGIIKIJJHJILILLLNLIGJLKKGMHMJGK@J9IIJGG+LIIGJJJE?8FFIIGK+DCJIDKB4GDFF1*9A@C@*0+	AS:i:-5	XS:i:-5	XN:i:0	XM:i:2	XO:i:0	XG:i:	NM:i:2	MD:Z:91C7C0	YS:i:0	YN:i:-60	Yn:i:0	ZN:i:-60	Zn:i:0	YT:Z:CP	Zm:Z:4	Zp:Z:5.088''')
                b = parse_sam_loc_mapq(
                    '''ERR050082.24501	83	1	1639728	41	100M	=	1639307	-521	GGTGCTGCCACAGGCAGGATGCGGGCTCCGCTTCAGCTAAGCAACAAGTGTTCCCAAGAATGGATATGGAGGCTGGGCGCGGTGGCTCACGCCTGTAATC	&=7@@9CEC?I:KGFDFBJ?GBJHIK@JKEHKJHJLKI?ELIGHKIJIAIKJFLKGJJEJGKLKEGMNMNKLNLLJJILCNKMMIKKIDJIKKJKHFHBD	NM:i:1	MD:Z:0A99	AS:i:99	XS:i:89	XA:Z:1,-1576518,100M,3;	Zm:Z:27	Zp:Z:40.752''')
                c = parse_sam_loc_mapq(
                    '''ERR050082.24501	83	1	1639728	41	100M	=	1639307	-521	GGTGCTGCCACAGGCAGGATGCGGGCTCCGCTTCAGCTAAGCAACAAGTGTTCCCAAGAATGGATATGGAGGCTGGGCGCGGTGGCTCACGCCTGTAATC	&=7@@9CEC?I:KGFDFBJ?GBJHIK@JKEHKJHJLKI?ELIGHKIJIAIKJFLKGJJEJGKLKEGMNMNKLNLLJJILCNKMMIKKIDJIKKJKHFHBD	NM:i:1	MD:Z:0A99	AS:i:99	XS:i:89	XA:Z:1,-1576518,100M,3;	Zm:Z:27	Zp:Z:40.752''')
                sim = same_alns_matrix([a, b, c], wiggle=10)
                self.assertTrue(numpy.all(numpy.array([[True,  False, False],
                                                       [False, True,  True],
                                                       [False, True,  True]]) == sim))
                sim = same_alns_matrix([a, b, c], wiggle=1000)
                self.assertTrue(numpy.all(numpy.array([[True, True, True],
                                                       [True, True, True],
                                                       [True, True, True]]) == sim))

            def test_same_alns_tiered_1(self):
                a = parse_sam_loc_mapq(
                    '''ERR050082.24499	163	1	1639199	5	100M	=	1639570	471	CGGCTGAGACAGAGCCCGGATGCTGAGCTGGGAGGAGGCGTCGGGTGTCATGTGGGGGACAAGCCCACATCCACGTCCACCAGGCTGAGGAAATAACCTA	:FDGEHFFCAKIJGCKMJJHJDGIIKIJJHJILILLLNLIGJLKKGMHMJGK@J9IIJGG+LIIGJJJE?8FFIIGK+DCJIDKB4GDFF1*9A@C@*0+	AS:i:-5	XS:i:-5	XN:i:0	XM:i:2	XO:i:0	XG:i:	NM:i:2	MD:Z:91C7C0	YS:i:0	YN:i:-60	Yn:i:0	ZN:i:-60	Zn:i:0	YT:Z:CP	Zm:Z:4	Zp:Z:5.088''')
                b = parse_sam_loc_mapq(
                    '''ERR050082.24501	83	1	1639728	41	100M	=	1639307	-521	GGTGCTGCCACAGGCAGGATGCGGGCTCCGCTTCAGCTAAGCAACAAGTGTTCCCAAGAATGGATATGGAGGCTGGGCGCGGTGGCTCACGCCTGTAATC	&=7@@9CEC?I:KGFDFBJ?GBJHIK@JKEHKJHJLKI?ELIGHKIJIAIKJFLKGJJEJGKLKEGMNMNKLNLLJJILCNKMMIKKIDJIKKJKHFHBD	NM:i:1	MD:Z:0A99	AS:i:99	XS:i:89	XA:Z:1,-1576518,100M,3;	Zm:Z:27	Zp:Z:40.752''')
                c = parse_sam_loc_mapq(
                    '''ERR050082.24501	83	1	1639728	41	100M	=	1639307	-521	GGTGCTGCCACAGGCAGGATGCGGGCTCCGCTTCAGCTAAGCAACAAGTGTTCCCAAGAATGGATATGGAGGCTGGGCGCGGTGGCTCACGCCTGTAATC	&=7@@9CEC?I:KGFDFBJ?GBJHIK@JKEHKJHJLKI?ELIGHKIJIAIKJFLKGJJEJGKLKEGMNMNKLNLLJJILCNKMMIKKIDJIKKJKHFHBD	NM:i:1	MD:Z:0A99	AS:i:99	XS:i:89	XA:Z:1,-1576518,100M,3;	Zm:Z:27	Zp:Z:40.752''')
                better_tiers = [[], [], []]
                equal_tiers = [[0, 1, 2], [0, 1, 2], [0, 1, 2]]
                sim = same_alns_tiered([a, b, c], better_tiers, equal_tiers, wiggle=10)
                self.assertEqual((0, 1, 1, True, False), sim[0])
                self.assertEqual((0, 2, 2, True, False), sim[1])
                self.assertEqual((0, 2, 2, True, False), sim[2])
                sim = same_alns_tiered([a, b, c], better_tiers, equal_tiers, wiggle=1000)
                self.assertEqual((0, 3, 3, True, True), sim[0])
                self.assertEqual((0, 3, 3, True, True), sim[1])
                self.assertEqual((0, 3, 3, True, True), sim[2])

            def test_same_alns_tiered_2(self):
                a = parse_sam_loc_mapq(
                    '''ERR050082.24499	83	1	1639728	5	100M	=	1639570	471	CGGCTGAGACAGAGCCCGGATGCTGAGCTGGGAGGAGGCGTCGGGTGTCATGTGGGGGACAAGCCCACATCCACGTCCACCAGGCTGAGGAAATAACCTA	:FDGEHFFCAKIJGCKMJJHJDGIIKIJJHJILILLLNLIGJLKKGMHMJGK@J9IIJGG+LIIGJJJE?8FFIIGK+DCJIDKB4GDFF1*9A@C@*0+	AS:i:-5	XS:i:-5	XN:i:0	XM:i:2	XO:i:0	XG:i:	NM:i:2	MD:Z:91C7C0	YS:i:0	YN:i:-60	Yn:i:0	ZN:i:-60	Zn:i:0	YT:Z:CP	Zm:Z:4	Zp:Z:5.088''')
                b = parse_sam_loc_mapq(
                    '''ERR050082.24501	83	1	1639728	41	100M	=	1639307	-521	GGTGCTGCCACAGGCAGGATGCGGGCTCCGCTTCAGCTAAGCAACAAGTGTTCCCAAGAATGGATATGGAGGCTGGGCGCGGTGGCTCACGCCTGTAATC	&=7@@9CEC?I:KGFDFBJ?GBJHIK@JKEHKJHJLKI?ELIGHKIJIAIKJFLKGJJEJGKLKEGMNMNKLNLLJJILCNKMMIKKIDJIKKJKHFHBD	NM:i:1	MD:Z:0A99	AS:i:99	XS:i:89	XA:Z:1,-1576518,100M,3;	Zm:Z:27	Zp:Z:40.752''')
                c = parse_sam_loc_mapq(
                    '''ERR050082.24501	83	1	1639728	41	100M	=	1639307	-521	GGTGCTGCCACAGGCAGGATGCGGGCTCCGCTTCAGCTAAGCAACAAGTGTTCCCAAGAATGGATATGGAGGCTGGGCGCGGTGGCTCACGCCTGTAATC	&=7@@9CEC?I:KGFDFBJ?GBJHIK@JKEHKJHJLKI?ELIGHKIJIAIKJFLKGJJEJGKLKEGMNMNKLNLLJJILCNKMMIKKIDJIKKJKHFHBD	NM:i:1	MD:Z:0A99	AS:i:99	XS:i:89	XA:Z:1,-1576518,100M,3;	Zm:Z:27	Zp:Z:40.752''')
                better_tiers = [[], [0], [0]]
                equal_tiers = [[0], [1, 2], [1, 2]]
                sim = same_alns_tiered([a, b, c], better_tiers, equal_tiers, wiggle=10)
                self.assertEqual((0, 1, 3, True, True), sim[0])
                self.assertEqual((1, 2, 3, True, True), sim[1])
                self.assertEqual((1, 2, 3, True, True), sim[2])

        unittest.main(argv=[sys.argv[0]])

    add_args(_parser)

    # Some basic flags
    _parser.add_argument('--profile', action='store_const', const=True, default=False, help='Print profiling info')
    _parser.add_argument('--verbose', action='store_const', const=True, default=False, help='Be talkative')
    _parser.add_argument('--version', action='store_const', const=True, default=False, help='Print version and quit')
    _parser.add_argument('--upto', metavar='int', type=int, help='Stop after this many alignments')

    _args = _parser.parse_args(sys.argv[1:])

    go_profile(vars(_args))
