import os
import re
import random
import sys
import string
from math import isnan
try:
    import numpypy as numpy
except ImportError:
    import numpy
from randutil import WeightedRandomGenerator
from read import Read
from reference import iter_fasta_chunks_fast
from collections import defaultdict

_revcomp_trans = string.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s):
    return s[::-1].translate(_revcomp_trans)


class FragmentSimSerial(object):
    """ Given a set of FASTA references, provides a generator that will
        generate a given number of fragment sequences, always from the
        forward strand, with a given fragment length.

        Downside is that it has to scan all of the FASTA, even when the
        user requests only a few reads.  Another downside is that it's
        difficult to generate a mix of read lengths over the course of
        one FASTA scan.

        Major advantage is that it uses very little memory. """

    def _estimate_length(self):
        """ Estimate the total length of the reference.  This is
            quicker than scanning all the FASTA files and usually
            provides an estimate that's too high by < 5%. """
        return sum([os.path.getsize(fn) for fn in self.fasta_filenames])

    def __init__(self, fasta_filenames):
        self.__re = re.compile('[^ACGTacgt]')  # detects ambiguous chars
        self.fasta_filenames = fasta_filenames
        self.total_length = self._estimate_length() * 0.9

    def simulate_batch(self, num_target, length, chunk_size=5000000):
        """ Simulate a batch of around 'num_target' fragments, each
            with length 'max_length'.  Use the given chunk size, which
            governs the amount of FASTA that gets parsed at a time. """
        last_seq, last_ref_id = '', None

        def my_binom(n_possible):
            return numpy.random.binomial(num_target, 1.0 * n_possible / self.total_length)

        for ref_id, offset, seq in iter_fasta_chunks_fast(self.fasta_filenames, chunk_size=chunk_size):
            if last_ref_id is None or ref_id != last_ref_id:
                last_ref_id = ref_id  # first chunk from this ref id
                ln_added = 0
            else:
                # grab a little bit from the previous chunk
                ln_added = min(length-1, len(last_seq))
                seq = last_seq[-ln_added:] + seq
                offset -= ln_added
            actual_length = len(seq) - (length - 1) + ln_added
            last_seq = seq
            if actual_length <= 0:
                continue
            tot_cov = my_binom(actual_length)
            if tot_cov == 0:
                continue
            # pick fragment starting positions
            for i in iter(numpy.random.random_integers(0, actual_length-1, tot_cov)):
                substr = seq[i:i+length]
                if len(substr) < length:
                    continue  # hit edge of ref_id
                if self.__re.search(substr):
                    continue  # ambiguous characters
                yield ref_id, offset + i, substr


class FragmentSimSerial2(object):

    def _estimate_length(self):
        """ Estimate the total length of the reference.  This is
            quicker than scanning all the FASTA files and usually
            provides an estimate that's too high by < 5%. """
        return sum([os.path.getsize(fn) for fn in self.fasta_filenames])

    def __init__(self,
                 fasta_filenames,  # FASTA fiels to simulate from
                 dists):  # distributions describing input data
        self.__re = re.compile('[^ACGTacgt]')  # detects ambiguous chars
        self.fasta_filenames = fasta_filenames
        self.dists = dists
        self.total_length = self._estimate_length() * 0.9
        self.conc_fraglen_avg = dists.sc_dist_conc.avg_fraglen
        self.disc_fraglen_avg = dists.sc_dist_disc.avg_fraglen
        self.unp_fraglen_avg = dists.sc_dist_unp.avg_fraglen
        self.bad_end_fraglen_avg = dists.sc_dist_bad_end.avg_fraglen
        self.max_fraglen = max(dists.longest_fragment(), dists.longest_unpaired())
        assert self.max_fraglen > 0

    def simulate_batch(self, fraction, minimum, chunk_size=5000000):
        """ Simulate a batch of around 'num_target' fragments.  Some of the
            fragments will be for concordant pairs, some for discordant pairs,
            and some for unpaired reads.  The mix is according to the
            frac_concordant and frac_discordant parameters given to the
            contsructor. """
        assert chunk_size >= self.max_fraglen
        num_concordant, num_discordant, num_unpaired, num_bad_end = 0, 0, 0, 0
        if not self.dists.sc_dist_conc.empty():
            num_added = self.dists.sc_dist_conc.num_added
            num_concordant = int(max(fraction * num_added, minimum))
        if not self.dists.sc_dist_disc.empty():
            num_added = self.dists.sc_dist_disc.num_added
            num_discordant = int(max(fraction * num_added, minimum))
        if not self.dists.sc_dist_unp.empty():
            num_added = self.dists.sc_dist_unp.num_added
            num_unpaired = int(max(fraction * num_added, minimum))
        if not self.dists.sc_dist_bad_end.empty():
            num_added = self.dists.sc_dist_bad_end.num_added
            num_bad_end = int(max(fraction * num_added, minimum))
        assert num_concordant + num_discordant + num_unpaired + num_bad_end > 0
        assert num_concordant >= 0 and num_discordant >= 0 and num_unpaired >= 0 and num_bad_end >= 0
        last_seq, last_ref_id = '', None

        def my_binom(length, has_spillover):
            """ Do three draws for the three types of fragments """
            if has_spillover:
                n_concordant_chances = n_discordant_chances = \
                    n_unpaired_chances = n_bad_end_chances = length - self.max_fraglen + 1
            else:
                n_concordant_chances = length - self.conc_fraglen_avg + 1
                n_bad_end_chances = length - self.bad_end_fraglen_avg + 1
                n_discordant_chances = length - self.disc_fraglen_avg + 1
                n_unpaired_chances = length - self.unp_fraglen_avg + 1
            return (0 if (num_concordant == 0 or n_concordant_chances <= 0) else numpy.random.binomial(num_concordant, float(n_concordant_chances) / self.total_length),
                    0 if (num_discordant == 0 or n_discordant_chances <= 0) else numpy.random.binomial(num_discordant, float(n_discordant_chances) / self.total_length),
                    0 if (num_unpaired == 0 or n_unpaired_chances <= 0) else numpy.random.binomial(num_unpaired, float(n_unpaired_chances) / self.total_length),
                    0 if (num_bad_end == 0 or n_bad_end_chances <= 0) else numpy.random.binomial(num_bad_end, float(n_bad_end_chances) / self.total_length))

        for ref_id, offset, seq in iter_fasta_chunks_fast(self.fasta_filenames, chunk_size=chunk_size):
            spillover = True
            if last_ref_id is None or ref_id != last_ref_id:
                last_ref_id = ref_id  # first chunk from this ref id
                spillover = False
            else:
                # add spillover from previous sequence
                assert len(last_seq) >= self.max_fraglen-1
                seq = last_seq[-(self.max_fraglen - 1):] + seq
                offset -= (self.max_fraglen - 1)
            #
            # ------------============================================
            #  spillover    new sequence
            #
            # len(spillover) == self.max_fraglen-1 or 0, can't be in between
            #
            # Say we're sampling a sequence with length self.max_fraglen - 1
            # There's one LESS spillover position it could come from, but one
            # MORE at the right.
            #
            # Now what about the no-spillover case.  Now say we're sampling a
            # sequence with length self.max_fraglen - 1.  There's on
            #
            nfrags_max = len(seq) - self.max_fraglen + 1
            last_seq = seq
            concordant_cov, discordant_cov, unpaired_cov, bad_end_cov = my_binom(len(seq), spillover)
            if concordant_cov == discordant_cov == unpaired_cov == bad_end_cov == 0:
                continue
            for cov, dist, typ in zip([concordant_cov, discordant_cov, unpaired_cov, bad_end_cov],
                                      [self.dists.sc_dist_conc, self.dists.sc_dist_disc,
                                       self.dists.sc_dist_unp, self.dists.sc_dist_bad_end],
                                      ['conc', 'disc', 'unp', 'bad_end']):
                if cov == 0:
                    continue
                unpaired = typ in ['unp', 'bad_end']
                # pick fragment starting positions
                draws = [dist.draw() for _ in xrange(cov)]
                # have to special-case this for unpaired draws
                fraglens = [int(round(d[5 if unpaired else 2])) for d in draws]
                assert not any(x > self.max_fraglen for x in fraglens)
                if spillover:
                    leftmost = [self.max_fraglen - fraglen for fraglen in fraglens]
                    rand_max = [nfrags_max] * cov
                else:
                    leftmost = [0] * cov
                    diffs = [self.max_fraglen - fraglen for fraglen in fraglens]
                    rand_max = [nfrags_max + diff for diff in diffs]
                for i, seq_off in enumerate(map(int, numpy.random.random(cov) * rand_max)):
                    left_off = leftmost[i] + seq_off
                    if fraglens[i] == 0:
                        fraglens[i] = self.max_fraglen
                    substr = seq[left_off:left_off+fraglens[i]]
                    if len(substr) < fraglens[i]:
                        continue  # hit edge of ref_id
                    if self.__re.search(substr):
                        continue  # ambiguous characters
                    ref_off = offset + left_off
                    if typ == 'unp':
                        # unpaired
                        sc_draw = draws[i]
                        fw, _, _, rf_aln, sc, rl, _, _ = sc_draw
                        assert rl == fraglens[i]
                        rdseq = substr if fw else revcomp(substr)
                        read = Read.from_simulator(rdseq, None, ref_id, ref_off, fw, sc, typ)
                        mutate(read, fw, sc_draw)  # mutate unpaired read
                        yield typ, read, None
                    elif typ == 'bad_end':
                        sc_draw = draws[i]
                        fw, _, _, rf_aln, sc, rl, mate1, ordlen = sc_draw
                        assert rl == fraglens[i]
                        rdseq = substr if fw else revcomp(substr)
                        rdp1 = Read.from_simulator(rdseq, None, ref_id, ref_off, fw, sc, typ)
                        mutate(rdp1, fw, sc_draw)  # mutate unpaired read
                        rdseq2 = ''.join([random.choice('ACGT') for _ in xrange(ordlen)])
                        # TODO: borrow qualities rather than make them up
                        qual2 = 'I' * len(rdseq2)
                        rdp2 = Read(rdp1.name + '2', rdseq2, qual2)
                        if not mate1:
                            rdp1, rdp2 = rdp2, rdp1
                        yield typ, rdp1, rdp2
                    else:
                        # paired-end
                        sc1_draw, sc2_draw, _, upstream1 = draws[i]
                        m1fw, _, _, rf_aln_1, sc1, rl1, _, _ = sc1_draw
                        m2fw, _, _, rf_aln_2, sc2, rl2, _, _ = sc2_draw
                        assert max(rl1, rl2) <= fraglens[i]
                        if upstream1:
                            seq1, seq2 = substr[:rl1], substr[-rl2:]
                        else:
                            seq2, seq1 = substr[:rl2], substr[-rl1:]
                        if not m1fw:
                            seq1 = revcomp(seq1)
                        if not m2fw:
                            seq2 = revcomp(seq2)
                        # Now we have the Watson offset for one mate or the other,
                        # depending on which mate is to the left w/r/t Watson.
                        if upstream1:
                            refoff1 = ref_off
                            refoff2 = ref_off + fraglens[i] - rl2
                        else:
                            refoff1 = ref_off + fraglens[i] - rl1
                            refoff2 = ref_off
                        rdp1 = Read.from_simulator(seq1, None, ref_id, refoff1, m1fw, sc1, typ)
                        rdp2 = Read.from_simulator(seq2, None, ref_id, refoff2, m2fw, sc2, typ)
                        mutate(rdp1, m1fw, sc1_draw)
                        mutate(rdp2, m2fw, sc2_draw)
                        yield typ, rdp1, rdp2


class SequenceSimulatorSerialStratifiedAlignability(object):

    def _estimate_length(self):
        return sum([os.path.getsize(fn) for fn in self.fasta_filenames])

    def __init__(self, fasta_filenames, alignability):
        self.__re = re.compile('[^ACGTacgt]')
        self.alignability = alignability
        self.fasta_filenames = fasta_filenames
        self.total_length = self._estimate_length() * 0.9

    def simulate_batch(self, num_target, max_length, chunk_size=10*1024):
        last_ref_id, ala_iter = None, None
        freq, num_strata = self.alignability.frequencies, self.alignability.max_strata
        per_stratum_target = num_target / num_strata
        num = 0
        num_in_strata = [0] * num_strata
        last_seq = ''
        ala, poss = [0] * chunk_size, []
        for i in xrange(0, num_strata):
            poss.append([0] * chunk_size)
        for ref_id, offset, seq in iter_fasta_chunks_fast(self.fasta_filenames, chunk_size=chunk_size):
            if last_ref_id is None or ref_id != last_ref_id:
                if not self.alignability.has_ref_id(ref_id):
                    continue  # skip refs for which we don't also have alignability
                self.alignability.load_ref_id(ref_id)
                ala_cursor = 0
                while self.alignability.records[ala_cursor][1] <= offset:
                    ala_cursor += 1
                last_ref_id = ref_id
                ln = max_length - 1
            else:
                ln = min(max_length-1, len(last_seq))
                seq = last_seq[-ln:] + seq
                offset -= ln
            actual_length = len(seq) - ln
            assert actual_length <= len(seq) - (max_length - 1)
            last_seq = seq

            # Populate stratum_hist
            stratum_hist = defaultdict(int)
            ref_cursor = offset
            n_left = actual_length
            ala_cursor_orig = ala_cursor
            while ala_cursor < len(self.alignability.records) and n_left > 0:
                st, en, stratum = self.alignability.records[ala_cursor]
                diff = min(en - ref_cursor, n_left)
                stratum_hist[stratum] += diff
                if diff == en - ref_cursor:
                    ala_cursor += 1
                ref_cursor += diff
                n_left -= diff
            assert n_left == 0 or ala_cursor == len(self.alignability.records)

            tot_covs = [numpy.random.binomial(per_stratum_target, 1.0 * stratum_hist[s] / freq[s]) for s in xrange(num_strata)]
            if not any(tot_covs):
                continue

            # Populate poss
            ref_cursor = offset
            n_left = actual_length
            ala_cursor = ala_cursor_orig
            stratum_cursor = [0] * num_strata
            while ala_cursor < len(self.alignability.records) and n_left > 0:
                st, en, stratum = self.alignability.records[ala_cursor]
                diff = min(en - ref_cursor, n_left)
                poss[stratum][stratum_cursor[stratum]:stratum_cursor[stratum] + diff] = xrange(ref_cursor - offset, ref_cursor - offset + diff)
                stratum_cursor[stratum] += diff
                if diff == en - ref_cursor:
                    ala_cursor += 1
                ref_cursor += diff
                n_left -= diff
            assert n_left == 0 or ala_cursor == len(self.alignability.records)

            for stratum, tot_cov in enumerate(tot_covs):
                if tot_cov == 0:
                    continue
                for pos_i in numpy.random.random_integers(0, stratum_hist[stratum]-1, tot_cov):
                    i = poss[stratum][pos_i]
                    substr = seq[i:i+max_length]
                    if len(substr) < max_length:
                        continue
                    if self.__re.search(substr):
                        continue
                    num += 1
                    num_in_strata[stratum] += 1
                    yield ref_id, offset + i, substr


class SequenceSimulator(object):
    """ Encapsulates a Reference and allows user to sample fragments of
        specified length.  Simulated fragments are not permitted to overlap a
        non-A/C/G/T character in the reference.  There is no awareness of the
        repeat content of the genome.
    """
    
    def __init__(self, ref):
        self.__re = re.compile('[^ACGT]')
        self.ref = ref
        lens = {}
        for name in ref.names():
            lens[name] = ref.length(name)
        self.rnd = WeightedRandomGenerator(lens)

    def _pick_ref(self, ln):
        """ Pick a reference sequence to simulate from in a weighted random
            fashion """
        attempts = 0
        while True:
            refi = self.rnd.next()
            assert self.ref.hasName(refi)
            if self.ref.length(refi) >= ln:
                return refi  # success
            attempts += 1
            if attempts > 5:
                raise RuntimeError('Tried %d times to find reference with length at least %d' % (attempts, ln))

    def _pick_substring(self, ref_id, length):
        """ Pick a substring from given reference sequence of given length
            that doesn't overlap any non-A/C/G/Ts """
        ref_offset = random.randint(0, self.ref.length(ref_id) - length)
        seq = self.ref.get(ref_id, ref_offset, length).upper()
        while self.__re.search(seq):
            ref_offset = random.randint(0, self.ref.length(ref_id) - length)
            seq = self.ref.get(ref_id, ref_offset, length).upper()
        return seq, ref_offset

    def sim(self, length):
        """ Simulate a read of length ln.  Substrings overlapping any
            non-A/C/G/T characters are rejected. """

        # Pick a reference sequence in a weighted random fashion
        ref_id = self._pick_ref(length)

        # Pick a substring that doesn't overlap any non-A/C/G/Ts
        seq, ref_offset = self._pick_substring(ref_id, length)

        # Possibly reverse complement
        fw = random.random() < 0.5
        if not fw:
            seq = revcomp(seq)
        assert not 'N' in seq

        return ref_id, ref_offset, fw, seq


def mutate(rd, rdfw, sc_dist_draw):
    """ Given a read that already has the appropriate length (i.e. equal to #
        characters on the reference side of the alignment), take the alignment
        information contained in the scDistDraw object and modify rd to
        contain the same pattern of edits.  Modifies rd in place. """
    fw, qual, rd_aln, rf_aln, sc, _, _, _ = sc_dist_draw
    if rdfw != fw:
        qual, rd_aln, rf_aln = qual[::-1], rd_aln[::-1], rf_aln[::-1]
    rd.qual = qual  # Install qualities
    # Walk along stacked alignment, making corresponding changes to
    # read sequence
    seq = []
    i, rdi, rfi = 0, 0, 0
    while i < len(rd_aln):
        if rd_aln[i] == rf_aln[i] and rf_aln[i] != 'N':
            seq.append(rd.seq[rfi])
            rfi += 1
            rdi += 1
        elif rd_aln[i] == '-':
            rfi += 1
        elif rf_aln[i] == '-':
            seq.append(rd_aln[i])
            rdi += 1
        elif rd_aln[i] != rf_aln[i] and (rd_aln[i] == 'N'):
            seq.append('N')
            rfi += 1
            rdi += 1
        elif rf_aln[i] == 'N':
            seq.append(random.choice('ACGT'))
            rfi += 1
            rdi += 1
        elif rd_aln[i] != rf_aln[i]:
            assert rfi < len(rd.seq)
            oldc = rd.seq[rfi].upper()
            cs = ['A', 'C', 'G', 'T']
            assert oldc in cs, "oldc was %s" % oldc
            cs.remove(oldc)
            newc = random.choice(cs)
            seq.append(newc)
            rfi += 1
            rdi += 1
        else:
            raise RuntimeError('Should not get here')
        i += 1
    rd.seq = ''.join(seq)
    assert len(rd.seq) == len(rd.qual), "%s\n%s\n%s\n%s" % (rd_aln, rf_aln, rd.seq, rd.qual)

if __name__ == "__main__":
    
    import argparse
    from shutil import rmtree
    from tempfile import mkdtemp
    from alignability import AlignabilityStratifiedIndexed

    parser = argparse.ArgumentParser()
    parser.add_argument('--test', action='store_const', const=True, default=False, help='Do unit tests')

    args = parser.parse_args()
    
    if args.test:
        import unittest
        
        random.seed(100)
        
        class TestMutate(unittest.TestCase):
            
            def test_1(self):
                # No edits
                rd = Read('r1', 'ACGT', '!!!!')
                mutate(rd, True, (True, 'IIII', 'TTTT', 'TTTT', 10))
                self.assertEqual(rd.seq, 'ACGT')
                self.assertEqual(rd.qual, 'IIII')
            
            def test_2(self):
                # 1 mismatch
                rd = Read('r2', 'ACGT', '!!!!')
                mutate(rd, True, (True, 'IIII', 'TTTT', 'TTTA', 10))
                self.assertNotEqual(rd.seq, 'ACGT')
                self.assertEqual(rd.qual, 'IIII')
            
            def test_3(self):
                # 1 read gap
                rd = Read('r3', 'A' * 12, '!' * 12)
                # rd: TTTTT--TTTTT
                # rf: TTTTTTTTTTTT
                mutate(rd, True, (True, 'I' * 10, 'TTTTT--TTTTT', 'T' * 12, 10))
                self.assertEqual(rd.seq, 'A' * 10)
                self.assertEqual(rd.qual, 'I' * 10)
            
            def test_4(self):
                # 1 reference gap
                rd = Read('r4', 'C' * 10, '#' * 10)
                # rd: TTTTTTTTTTTT
                # rf: TTTTT--TTTTT
                mutate(rd, True, (True, 'J' * 12, 'TTTTTTTTTTTT', 'TTTTT--TTTTT', 10))
                self.assertEqual(rd.seq, 'CCCCCTTCCCCC')
                self.assertEqual(rd.qual, 'J' * 12)
            
            def test_5(self):
                # N in the reference
                rd = Read('r5', 'CCCCC', '#####')
                # rd: TTTTT
                # rf: TTNTT
                mutate(rd, True, (True, 'KKKKK', 'TTTTT', 'TTNTT', 10))
                self.assertEqual(rd.seq[:2], 'CC')
                self.assertEqual(rd.seq[-2:], 'CC')
                self.assertNotEqual(rd.seq[2], 'N')
                self.assertEqual(rd.qual, 'KKKKK')
            
            def test_6(self):
                # N in the read
                rd = Read('r6', 'CCCCC', '#####')
                # rd: TTNTT
                # rf: TTTTT
                mutate(rd, True, (True, 'KKKKK', 'TTNTT', 'TTTTT', 10))
                self.assertEqual(rd.seq, 'CCNCC')
                self.assertEqual(rd.qual, 'KKKKK')
        
        class TestSequenceSimulatorSerial(unittest.TestCase):

            def setUp(self):
                self.tmpdir = mkdtemp()
                self.fa_fn_1 = os.path.join(self.tmpdir, 'tmp1.fa')
                with open(self.fa_fn_1, 'w') as fh:
                    fh.write('''>short_name1 with some stuff after whitespace
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>short_name2 with some stuff after whitespace
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
''')

            def tearDown(self):
                rmtree(self.tmpdir)

            def test_1(self):
                sss = FragmentSimSerial([self.fa_fn_1])
                simreads = [x for x in sss.simulate_batch(100, 120, chunk_size=100)]
                self.assertGreater(len(simreads), 50)
                tot = defaultdict(int)
                for ref_id, _, _ in simreads:
                    tot[ref_id] += 1
                self.assertGreater(tot['short_name1'], tot['short_name2'])

        class TestSequenceSimulatorSerialStratifiedAlignability(unittest.TestCase):

            def setUp(self):
                self.tmpdir = mkdtemp()
                self.fa_fn_1 = os.path.join(self.tmpdir, 'tmp1.fa')
                with open(self.fa_fn_1, 'w') as fh:
                    fh.write('''>short_name1 with some stuff after whitespace
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>short_name2 with some stuff after whitespace
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
>short_name3 with some stuff after whitespace
ACACACACACACACACACACACACACACACACACACACACACACACACACACACACACAC
ACACACACACACACACACACACACACACACACACACACACACACACACACACACACACAC
>short_name4 with some stuff after whitespace
ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT
ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT
>short_name5 with some stuff after whitespace
AGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG
AGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG
>short_name6 with some stuff after whitespace
CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG
CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG
''')
                self.wig_fn_1 = os.path.join(self.tmpdir, 'tmp1.wig')
                with open(self.wig_fn_1, 'w') as fh:
                    # first record intentionall does not extend to the
                    # beginning or end of short_name1
                    fh.write('''# comment
short_name1 100   1100   1
short_name2 0   300   0.5
short_name3 0   120   0.3
short_name4 0   120   0.2
short_name5 0   120   0.08
short_name6 0   110   0.03
''')

            def tearDown(self):
                rmtree(self.tmpdir)

            def test_1(self):
                ala = AlignabilityStratifiedIndexed([self.wig_fn_1], max_strata=6)
                sssa = SequenceSimulatorSerialStratifiedAlignability([self.fa_fn_1], ala)
                simreads = [x for x in sssa.simulate_batch(3000, 50, chunk_size=100)]
                self.assertGreater(len(simreads), 2000)
                tot = defaultdict(int)
                for ref_id, _, _ in simreads:
                    tot[ref_id] += 1
                for d in xrange(1, 7):
                    self.assertTrue(400 < tot['short_name%d' % d] < 600)

        unittest.main(argv=[sys.argv[0]])
