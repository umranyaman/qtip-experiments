import sys
import re

__author__ = 'langmead'

"""
Given SAM files from multiple different read alignment tools, iterate through
each, assuming alignments are in the same order in each, and if all the tools
align the read, ask how many.

Imagine 3 tools and 3 separate boxplots.  For a given tool, the boxplots show
(a) the distribution of MAPQs for alignments where this tool agrees with the
other 2, (b) the same for when the tool agrees with only 1 of the other 2, (c)
the same for when the tool disagrees with both.

How would I quibble with plots like these?  Maybe by saying that disagreement
is a poor proxy for incorrectness.  Could we build confidence by also doing a
simulation?

Could we build confidence by taking a given tool, running it with strictly
more stringent sensitivity settings and seeing a more pronounced mean shift
for the predicted MAPQs than for the
"""


VERSION = '0.0.0'


def add_args(parser):
    # Overall arguments
    parser.add_argument('--sam', metavar='path', type=str, nargs='+', required=True,
                        help='SAM file(s) containing aligned reads')

mapq_re = re.compile('Zm:[iZ]:([0-9]+)')
mapq_prec_re = re.compile('Zp:Z:([.0-9]+)')
digre = re.compile('([01-9]+)')


def cigar_to_list(cigar_string):
    spl = digre.split(cigar_string)
    return zip(spl[2::2], map(int, spl[1::2]))


def parse_sam_loc_mapq(st):
    toks = st.split('\t')
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


def count_same_alns(ls, wiggle=50):
    sim = [1] * len(ls)
    for i in range(len(ls)):
        for j in range(i+1, len(ls)):
            if compare_2_alns(ls[i], ls[j], wiggle=wiggle):
                sim[i] += 1
                sim[j] += 1
    return sim


def go(args):
    sams = map(lambda x: open(x, 'rb'), args['sam'])
    while True:
        lns = []
        for sam in sams:
            while True:
                ln = sam.readline()
                if len(ln) == 0:
                    lns.append(None)
                    break
                elif ln[0] == '@':
                    continue
                lns.append(ln)
                break

        assert len(lns) == len(sams)
        if all(map(lambda x: x is None, lns)):
            break
        assert not any((lambda x: x is None, lns))
        parsed_sam = map(parse_sam_loc_mapq, lns)
        sim = count_same_alns(parsed_sam, args['wiggle'])
        i = 1
        for ps, si in zip(parsed_sam, sim):
            rname, pos, mapq, mapq_orig, mapq_prec = ps
            print (','.join(map(str, [i, si, mapq, mapq_orig, mapq_prec])))
            i += 1


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
        print('Tandem simulator, version ' + VERSION)
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
                sim = count_same_alns([a, b, c], wiggle=10)
                self.assertEqual([1, 2, 2], sim)
                sim = count_same_alns([a, b, c], wiggle=1000)
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
                sim = count_same_alns([a, b, c, d, e], wiggle=10)
                self.assertEqual([1, 2, 2, 2, 2], sim)
                sim = count_same_alns([a, b, c, d, e], wiggle=1000)
                self.assertEqual([3, 3, 3, 2, 2], sim)
                sim = count_same_alns([a, b, c, d, e], wiggle=20000)
                self.assertEqual([5, 5, 5, 5, 5], sim)

        unittest.main(argv=[sys.argv[0]])

    add_args(_parser)

    # Some basic flags
    _parser.add_argument('--profile', action='store_const', const=True, default=False, help='Print profiling info')
    _parser.add_argument('--verbose', action='store_const', const=True, default=False, help='Be talkative')
    _parser.add_argument('--version', action='store_const', const=True, default=False, help='Print version and quit')

    _args = _parser.parse_args(sys.argv[1:])

    go_profile(vars(_args))