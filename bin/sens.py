#!/usr/bin/env python

"""
sens.py

Given a tool, tool parameters, and a genome index, assess the tool's ability to
find best-stratum alignments when aligning reads to that genome.  The tool uses
pseudo-random numbers to vary the number and configuration of the mismatches,
and the number, configuration and lengths of the gaps.

(Do we try to infer the scoring scheme?)
(Do we introduce Ns?) 
"""

import sys
import os
import math
import random
import bisect
import string
import argparse
import global_align

parser = argparse.ArgumentParser(\
    description='Evaluate the sensitivity of an alignment tool w/r/t given'
                'genome and set of alignment parameters.')

parser.add_argument(\
    '--fasta', metavar='path', dest='fasta', type=str, nargs='+',
    required=True, help='FASTA file(s) containing reference genome sequences')
parser.add_argument(\
    '--scoring', metavar='str', type=str, required=False,
    default='1,2,6,1,5,3,5,3',
    help='MatchBonus,MismatchMinPen,MismatchMaxPen,NPen,ReadGapConst,ReadGapLinear,RefGapConst,RefGapLinear')
parser.add_argument(\
    '--num-reads', metavar='int', dest='num_reads', type=int, default=100,
    required=False, help='Number of reads to simulate')
parser.add_argument(\
    '--min-id', metavar='fraction', dest='minid', type=float, default=0.95,
    required=False, help='Minimum pct identity to allow')
parser.add_argument(\
    '--min-len', metavar='int', dest='minlen', type=int, default=50,
    required=False, help='Minimum read length to simulate')
parser.add_argument(\
    '--max-len', metavar='int', dest='maxlen', type=int, default=400,
    required=False, help='Maximum read length to simulate')
parser.add_argument(\
    '--sampling', metavar='scheme', type=str, default="score",
    required=False, help='"score" or "nedit"')
parser.add_argument(\
    '--sanity', dest='sanity', action='store_const', const=True, default=False,
    help='Do various sanity checks')
parser.add_argument(\
    '--test', dest='test', action='store_const', const=True, default=False,
    help='Do unit tests')
parser.add_argument(\
    '--profile', dest='profile', action='store_const', const=True,
    default=False, help='Print profiling info')
parser.add_argument(\
    '--verbose', dest='verbose', action='store_const', const=True,
    default=False, help='Be talkative')
parser.add_argument(\
    '--version', dest='version', action='store_const', const=True,
    default=False, help='Print version and quit')

args = parser.parse_args()

_revcomp_trans = string.maketrans("ACGT", "TGCA")
def revcomp(s):
    return s[::-1].translate(_revcomp_trans)

class Scoring(object):
    
    """ Parameters governing how to score differences in an alignment """ 
    
    def __init__(self, ma, mmp, np, rdg, rfg):
        self.ma = ma             # match bonus: integer
        self.mmp = mmp           # mismatch penalty: (min, max) integer tuple
        self.np = np             # N penalty: integer
        self.rdg = rdg           # affine read gap penalty: (a, b)
        self.rfg = rfg           # affine reference gap penalty: (a, b)
    
    def rfgPenalty(self, length=1):
        return self.rfg[0] + length * self.rfg[1]
    
    def rdgPenalty(self, length=1):
        return self.rdg[0] + length * self.rdg[1]

class SimpleFunc(object):
    
    """ One-argument function with constant term and term that is proportional
        to some function of the argument.  Function can be f(x) = 0 (constant),
        f(x) = x (linear), f(x) = ln(x) (log), or f(x) = sqrt(x) (sqrt) """
    
    def __init__(self, type="const", I=None, X=None, C=0.0, L=0.0):
        self.type = type
        self.I = I
        self.X = X
        self.C = C
        self.L = L
        if I is None: I = float('-inf')
        if X is None: X = float('inf')
    
    def f(self, x):
        if self.type == "const":
            x = 0.0
        elif self.type == "sqrt":
            x = math.sqrt(x)
        elif self.type == "log":
            x = math.log(x)
        elif self.type != "linear":
            raise RuntimeError("Bad simple function type: '%s'" % self.type)
        return min(X, max(I, self.C + self.L * x))

class WeightedRandomGenerator(object):
    
    """ Given an ordered list of weights, generate with each call to next() an
        offset into the list of the weights with probability equal to the
        fraction of the total weight. """

    def __init__(self, weights):
        self.totals = []
        running_total = 0
        for w in iter(weights):
            running_total += w
            self.totals.append(running_total)

    def next(self):
        rnd = random.random() * self.totals[-1]
        return bisect.bisect_right(self.totals, rnd)

class ReadLengthGen(object):
    
    """ Given a maximum and minimum read length, generate with each call to
        next() an integer read length sampled uniformly from that range """ 
    
    def __init__(self, mn=50, mx=400):
        self.mn = mn
        self.mx = mx
    
    def next(self):
        return random.randint(self.mn, self.mx)

class MyMutableString(object):
    
    """ A string vaguely like Python's built-in string type but where it's
        easier and more efficient to insert things in the middle and delete
        stretches. """
    
    def __init__(self, s=None):
        self.slist = []
        if s is not None:
            self.append(s)

    def __str__(self):
        return ''.join(self.slist)
    
    def append(self, s):
        self.slist.extend(list(s))
    
    def delete(self, i, j=None):
        if j is None:
            j = i + 1
        for k in xrange(i, j):
            self.slist[k] = ""
    
    def set(self, i, s):
        assert i < len(self.slist)
        self.slist[i] = s
    
    def get(self, i):
        return self.slist[i]

class Mutator(object):
    
    """ Class that, given a read sequence, returns a version of the read
        sequence mutated somehow, along with a distance (score) indicating how
        much it was mutated.
        
        Doesn't take quality values into account. """
    
    def __init__(self, scoring, sampling="score", minId=0.7, maxRdLen=400, affine=False, verbose=False):
        assert minId >= 0.0 and minId <= 1.0
        self.minId = minId
        self.maxRdLen = maxRdLen
        self.scoring = scoring
        self.sampling = sampling
        self.sc_min = float('inf')
        self.sc_max = float('-inf')
        self.scoreToCombo = dict()
        self.combos_by_nedit = []
        self.combos_by_score = []
        self.affine = affine
        self.setupCombos(verbose=verbose)
    
    def addCombo(self, combo):
        
        """ Add an association between the given combination of mismatches and
            gaps, and the score that the combination yields. """
        
        nmm, nrdg, nrfg = combo
        nedit = sum(combo)
        sc = nmm * (self.scoring.mmp[1]) + \
             nrdg * (self.scoring.rdg[0] + self.scoring.rdg[1]) + \
             nrfg * (self.scoring.rfg[0] + self.scoring.rfg[1])
        self.sc_min = min(self.sc_min, sc)
        self.sc_max = max(self.sc_max, sc)
        if sc not in self.scoreToCombo:
            self.scoreToCombo[sc] = []
        self.scoreToCombo[sc].append(combo)
        self.combos_by_nedit.append((nedit, sc, nmm, nrdg, nrfg))
        self.combos_by_score.append((sc, nedit, nmm, nrdg, nrfg))
        self.combos_min_score = []
        self.combos_max_score = []
    
    def setupCombos(self, verbose=False):
        
        """ Sweep through all possible combinations of # mismatches, # read
            gaps, read gap lengths, # ref gaps, ref gap length that come in
            under the specified % identity minimum. """
        
        max_edits = int(math.ceil(self.maxRdLen * (1.0 - self.minId)))
        max_mms = max_edits
        max_rdgaps = max_edits
        max_rfgaps = max_edits
        max_rdgaplen = max_edits
        max_rfgaplen = max_edits
        id = 1.0
        for i in xrange(0, max_mms+1):
            id1 = id - (float(i) / self.maxRdLen)
            if id1 < self.minId:
                break # id1 is only going to decrease as we go further
            sc_mm = self.scoring.mmp * i
            for j in xrange(0, max_rdgaps+1):
                id2 = id1 - (float(j) / self.maxRdLen)
                if id2 < self.minId:
                    break # id2 is only going to decrease as we go further
                sc_rdgap = self.scoring.rdg[0] + self.scoring.rdg[1] * j
                for k in xrange(0, max_rfgaps+1):
                    id3 = id2 - (float(k) / self.maxRdLen)
                    if id3 < self.minId:
                        break # id3 is only going to decrease as we go further
                    sc_rfgap = self.scoring.rfg[0] + self.scoring.rfg[1] * k
                    # Score what we have so far
                    self.addCombo((i, j, k))
        # Sort combos by # edits, then by score
        if verbose:
            print >>sys.stderr, "Sorting combos by # edits..."
        self.combos_by_nedit.sort()
        print >>sys.stderr, "...done"
        if verbose:
            print >>sys.stderr, "Sorting combos by score..."
        self.combos_by_score.sort()
        print >>sys.stderr, "...done"
        # Build a list parallel to self.combos_by_score containing maximum
        # score seen so far up to each point
        min_sc, max_sc = float('inf'), float('-inf')
        for i in xrange(0, len(self.combos_by_score)):
            min_sc = min(min_sc, self.combos_by_score[i][0])
            max_sc = max(max_sc, self.combos_by_score[i][0])
            print(min_sc, max_sc)
            self.combos_min_score.append(min_sc) 
            self.combos_max_score.append(max_sc) 
    
    def __addEdits(self, orig_rd, sc, nmm, nrdg, nrfg, maxAttempts=10, verbose=False):
        
        """ Given a substring extracted from the reference genome (s) and given
            a # of mismatches, # of read gaps, and # of reference gaps, mutate
            the string with the given number of gaps and mismatches. """ 
        
        # There are several subtleties here.  First, applying a bunch of edits
        # whose penalties add to the desired score does not necessarily result
        # in a read that aligns with the desired score.  One reason is that
        # multiple edits near each other might "cancel each other out", e.g.
        # adjacent read and reference gaps could be more optimally scored as a
        # single mismatch (or match!).  Or two adjacent read gaps can be more
        # optimally scored as a single gap open with several extensions.
        #
        # We choose a very simple (naive?) method for resolving this.  We
        # simply insert the edits in random places, then check whether we ended
        # up with a proper score by calling a global aligner.
        attempts = 0
        rd = None
        rf = orig_rd[:]
        while True:
            rd = MyMutableString(orig_rd)
            assert str(rd) == orig_rd
            attempts += 1
            if attempts > maxAttempts:
                return orig_rd, False
            if attempts > 1:
                print >> sys.stderr, "Warning: attempt #%d" % attempts
            # Generate mismatches and gaps and place them in random spots along
            # the read
            mmset, rdgset, rfgset = set(), set(), set()
            while len(mmset) < nmm:
                # Insert a mismatch
                mmset.add(random.randint(0, len(orig_rd)-1))
            while len(rdgset) < nrdg:
                # Insert a read gap
                rdgset.add(random.randint(0, len(orig_rd)-1))
            while len(rfgset) < nrfg:
                # Insert a reference gap (prior to the character at the given
                # index)
                rfgset.add(random.randint(1, len(orig_rd)-1))
            
            # Apply mismatches to the read string
            for mm in mmset:
                origc = rd.get(mm)
                assert len(origc) == 1
                c = random.choice(['A', 'C', 'G', 'T'])
                while c == origc:
                    c = random.choice(['A', 'C', 'G', 'T'])
                rd.set(mm, c)
            
            # Apply read gaps, which involve inserting characters into the
            # reference string
            for gapi in sorted(list(rdgset) + list(rfgset), reverse=True):
                if gapi in rdgset:
                    rd.delete(gapi)
                if gapi in rfgset:
                    c = random.choice(['A', 'C', 'G', 'T'])
                    rd.set(gapi, c + rd.get(gapi))
            
            rdstr = str(rd)
            if not self.affine:
                scg, _, _ = global_align.globalAlign(rdstr, rf, lambda x, y: self.scoring.mmp[1] if x != y else 0, self.scoring.rdg[0] + self.scoring.rdg[1], self.scoring.rfg[0] + self.scoring.rfg[1], backtrace=False)
            else:
                scg, _, _ = global_align.globalAlignAffine(rdstr, rf, lambda x, y: self.scoring.mmp[1] if x != y else 0, self.scoring.rdg[0], self.scoring.rdg[1], self.scoring.rfg[0], self.scoring.rfg[1], backtrace=False)

            print "Shooting for %d, got %d (%d mismatches, %d read gaps, %d ref gaps)" % (sc, scg, nmm, nrdg, nrfg)
            if scg == sc:
                # Success!
                break
        # Check w/ global alignment that the final answer has expected score
        return rd, True
    
    def mutate(self, rd, verbose=False):
        
        """ First we find some (perhaps all) combinations that don't fall below
            the minimum identity threshold.  From among these, we pick one
            uniformly at random and try to apply it to the read.  This may or
            may not succeed (depending, e.g. on where gaps land). """
        
        if verbose:
            print >>sys.stderr, "Mutate called..."
        max_edits = int(math.ceil(len(rd) * (1.0 - self.minId)))
        max_elt = bisect.bisect_left(self.combos_by_nedit, (max_edits+1, None, None, None, None))
        assert max_elt > 0
        if self.sampling == "nedit": 
            assert max_elt > 0
            # Pick from among combos with acceptable number of edits,
            # uniformly at random.
            combo = self.combos_by_nedit[random.randint(0, max_elt-1)]
            (nedit, sc, nmm, nrdg, nrfg) = combo
            print ("selected combo: " + str(combo))
            # TODO: Actually mutate the read sequence
            return (rd, sc)
        elif self.sampling == "score":
            min_sc, max_sc = self.combos_min_score[max_elt-1], self.combos_max_score[max_elt-1]
            attempt = 0
            mutrd, sc, succ = None, None, None
            while not succ:
                attempt += 1
                sc = random.randint(min_sc, max_sc)
                min_sc_i = bisect.bisect_left(self.combos_by_score, (sc, None, None, None, None))
                min_sc_f = bisect.bisect_left(self.combos_by_score, (sc+1, None, None, None, None))
                assert min_sc_i <= min_sc_f
                if min_sc_f == min_sc_i:
                    # No scores actually started with sc, pick new sc to be next
                    # highest score that does exist in the list
                    sc = self.combos_by_score[min_sc_i][0]
                    min_sc_i = bisect.bisect_left(self.combos_by_score, (sc, None, None, None, None))
                    min_sc_f = bisect.bisect_left(self.combos_by_score, (sc+1, None, None, None, None))
                assert min_sc_i < len(self.combos_by_score)
                assert min_sc_f > min_sc_i, "sc=%d, from [minsc=%d, maxsc=%d]; min_sc top, bot=%d, %d, len(combos)=%d" % (sc, min_sc, max_sc, min_sc_i, min_sc_f, len(self.combos_by_score))
                rndi = random.randint(min_sc_i, min_sc_f-1)
                print "Combo: %s, read len: %d, attempt: %d" % (str(self.combos_by_score[rndi]), len(rd), attempt)
                sc, nedit, nmm, nrdg, nrfg = self.combos_by_score[rndi]
                mutrd, succ = self.__addEdits(rd, sc, nmm, nrdg, nrfg)
            assert sc is not None
            assert mutrd is not None
            return (mutrd, sc)
        else:
            raise RuntimeError("Unknown sampling scheme '%s'" % self.sampling)

class Simulator(object):
    
    """ Class that, given a collection of FASTA files, samples intervals of
        specified length from the strings contained in them. """
    
    def __init__(self, fafns):
        self.refs = dict()
        self.names = []
        self.lens = []
        for fafn in fafns:
            fafh = open(fafn, 'r')
            name = None
            for line in fafh:
                line = line.rstrip()
                if line.startswith('>'):
                    ind = line.index(" ")
                    if ind == -1: ind = len(line)
                    line = line[1:ind]
                    name = line
                    self.refs[name] = ""
                    self.names.append(name)
                    self.lens.append(0)
                else:
                    assert name is not None
                    self.refs[name] += line
                    self.lens[-1] += len(line)
            fafh.close()
        self.rnd = WeightedRandomGenerator(self.lens)
    
    def sim(self, ln, verbose=False):
        if verbose:
            print >>sys.stderr, "sim called..."
        refi = self.rnd.next()
        assert refi < len(self.names)
        off = random.randint(0, self.lens[refi] - ln)
        fw = True
        nm = self.names[refi]
        seq = self.refs[nm][off:off+ln]
        if random.random() > 0.5:
            fw = False
            seq = revcomp(seq)
        if verbose:
            print >>sys.stderr, "...done"
        return (nm, off, fw, seq)

class CaseGen(object):
    
    def __init__(self, fafns, scoring, sampling, minlen, maxlen, minId, verbose=False):
        self.minlen = minlen
        self.maxlen = maxlen
        self.scoring = scoring
        self.verbose = verbose
        self.rlg = ReadLengthGen(mn=minlen, mx=maxlen)
        self.sim = Simulator(args.fasta)
        self.mut = Mutator(scoring, sampling, minId, maxlen, verbose)

    def __iter__(self):
        return self
    
    def next(self):
        ln = self.rlg.next()
        ref, off, fw, seq = self.sim.sim(ln, verbose=self.verbose)
        mseq, sc = self.mut.mutate(seq, verbose=self.verbose)
        return seq

def go():
    sctok = map(int, string.split(args.scoring, ','))
    scoring = Scoring(sctok[0], (sctok[1], sctok[2]), sctok[3], (sctok[4], sctok[5]), (sctok[6], sctok[7]))
    cg = CaseGen(args.fasta, scoring, args.sampling, args.minlen, args.maxlen, args.minid, verbose=args.verbose)
    for x in xrange(0, args.num_reads):
        c = cg.next()

if args.profile:
    import cProfile
    cProfile.run('go()')
else:
    go()
