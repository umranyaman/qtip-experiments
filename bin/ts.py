#!/usr/bin/env python

"""
ts.py

A "tandem simulator," which wraps an alignment tool as it runs, eavesdrops on
the input and output, and builds a model that can be used to improve the
quality values calculated for aligned reads.

Things we learn from reads
==========================

- Read length distribution
- Quality values

Things we learn from alignments
===============================

- Alignment type (aligned, unaligned, concordant, discordant)
- Fragment length distribution
- Number and placement of mismatches and gaps

"""

__author__ = "Ben Langmead"
__email__ = "langmea@cs.jhu.edu"

import os
import sys
import re
import threading
import math
import string
import random
import bisect
import subprocess
import tempfile
import signal
import traceback
import cPickle
import time
from Queue import Queue
from sklearn.neighbors import KNeighborsClassifier

def quit_handler(signum,frame):
    traceback.print_stack()

signal.signal(signal.SIGQUIT,quit_handler)

_revcomp_trans = string.maketrans("ACGTacgt", "TGCAtgca")
def revcomp(s):
    return s[::-1].translate(_revcomp_trans)

class Read(object):
    """ Encapsulates one read """
    def __init__(self, name, seq, qual, orig=None):
        self.name = name
        self.seq = seq
        self.qual = qual
        self.orig = orig
        assert self.repOk()
    
    @classmethod
    def fromSimulator(cls, seq, qual, refid, refoff, fw, sc):
        # Construct appropriate name
        rdname = "!!ts-sep!!".join(["!!ts!!", refid, "+" if fw else "-", str(refoff), str(sc)])
        return cls(rdname, seq, qual)
    
    def __len__(self):
        """ Return number of nucleotides in read """
        return len(self.seq)
    
    def __str__(self):
        """ Return string representation """
        if self.orig is not None:
            return self.orig # original string preferred
        elif self.qual is not None:
            assert len(self.seq) == len(self.qual)
            return "@%s\n%s\n+\n%s\n" % (self.name, self.seq, self.qual)
        else:
            return ">%s\n%s\n" % (self.name, self.seq)
    
    def repOk(self):
        """ Check that read is internally consistent """
        if self.qual is not None:
            assert len(self.seq) == len(self.qual)
        return True

class Alignment(object):
    """ Encapsulates an alignment record for a single aligned read """

    __asRe  = re.compile('AS:i:([-]?[0-9]+)')
    __xsRe  = re.compile('XS:i:([-]?[0-9]+)')
    __ytRe  = re.compile('YT:Z:([A-Z]+)')
    __ysRe  = re.compile('YS:i:([-]?[0-9]+)')
    __mdRe  = re.compile('MD:Z:([^\s]+)')
    __xlsRe = re.compile('Xs:i:([-]?[0-9]+)')
    
    def __init__(self, ln):
        self.name, self.flags, self.refid, self.pos, self.mapq, self.cigar, \
        self.rnext, self.pnext, self.tlen, self.seq, self.qual, self.extra = \
        string.split(ln, '\t', 11)
        assert self.flags != "*"
        assert self.pos != "*"
        assert self.mapq != "*"
        assert self.tlen != "*"
        self.flags = int(self.flags)
        self.pos = int(self.pos)
        self.mapq = int(self.mapq)
        self.tlen = int(self.tlen)
        self.fw = (self.flags & 16) == 0
        self.mate1 = (self.flags & 64) != 0
        self.mate2 = (self.flags & 128) != 0
        self.paired = self.mate1 or self.mate2
        se = Alignment.__asRe.search(self.extra)
        self.bestScore = None
        if se is not None:
            self.bestScore = int(se.group(1))
        se = Alignment.__xsRe.search(self.extra)
        self.secondBestScore = None
        if se is not None:
            self.secondBestScore = int(se.group(1))
        se = Alignment.__xlsRe.search(self.extra)
        self.thirdBestScore = None
        if se is not None:
            self.thirdBestScore = int(se.group(1))
        se = Alignment.__ytRe.search(self.extra)
        self.alType = None
        if se is not None:
            self.alType = se.group(1)
        se = Alignment.__ysRe.search(self.extra)
        self.mateBest = None
        if se is not None:
            self.mateBest = int(se.group(1))
        self.mdz = None
        se = Alignment.__mdRe.search(self.extra)
        if se is not None:
            self.mdz = se.group(1)
        self.concordant = self.alType == "CP"
        self.discordant = self.alType == "DP"
        self.unpPair = self.alType == "UP"
        self.unp = self.alType == "UU"
        assert self.repOk()
    
    def isAligned(self):
        """ Return true iff read aligned """
        return (self.flags & 4) == 0
    
    def orientation(self):
        """ Return orientation as + or - """
        if (self.flags & 16) != 0:
            return "-"
        else:
            return "+"
    
    def fragmentLength(self):
        """ Return fragment length """
        return abs(self.tlen)
    
    def __len__(self):
        """ Return read length """
        return len(self.seq)
    
    def repOk(self):
        """ Check alignment for internal consistency """
        assert self.alType is not None
        assert self.paired or self.fragmentLength() == 0
        assert not self.isAligned() or self.bestScore is not None
        assert self.alType in ["CP", "DP", "UP", "UU"]
        return True

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

def cigarToList(cigar):
    """ Parse CIGAR string into a list of CIGAR operations. """
    ret = [];
    i = 0
    # CIGAR operations: MIDNSHP
    op_map = {'M':0, '=':0, 'X':0, 'I':1, 'D':2, 'N':3, 'S':4, 'H':5, 'P':6}
    while i < len(cigar):
        run = 0
        while i < len(cigar) and cigar[i].isdigit():
            run *= 10
            run += int(cigar[i])
            i += 1
        assert i < len(cigar)
        op = cigar[i]
        i += 1
        assert op in op_map
        ret.append([op_map[op], run])
    return ret

def mdzToList(md):
    """ Parse MD:Z string into a list of operations, where 0=match,
        1=read gap, 2=mismatch. """
    i = 0;
    ret = [] # list of (op, run, str) tuples
    while i < len(md):
        if md[i].isdigit(): # stretch of matches
            run = 0
            while i < len(md) and md[i].isdigit():
                run *= 10
                run += int(md[i])
                i += 1 # skip over digit
            if run > 0:
                ret.append([0, run, ""])
        elif md[i].isalpha(): # stretch of mismatches
            mmstr = ""
            while i < len(md) and md[i].isalpha():
                mmstr += md[i]
                i += 1
            assert len(mmstr) > 0
            ret.append([1, len(mmstr), mmstr])
        elif md[i] == "^": # read gap
            i += 1 # skip over ^
            refstr = ""
            while i < len(md) and md[i].isalpha():
                refstr += md[i]
                i += 1 # skip over inserted character
            assert len(refstr) > 0
            ret.append([2, len(refstr), refstr])
        else: assert False
    return ret

def cigarMdzToStacked(seq, cgp, mdp_orig):
    """ Takes parsed cigar and parsed MD:Z and generates a stacked alignment:
        a pair of strings with gap characters inserted (possibly) and where
        characters at at the same offsets are opposite each other in the
        alignment.  Returns tuple of 2 parallel strings: read string, ref
        string. """
    mdp = mdp_orig[:]
    rds, rfs = "", ""
    mdo, rdoff = 0, 0
    for c in cgp:
        assert mdo < len(mdp), "mdo=%d" % mdo
        op, run = c
        if op == 0:   # M
            # Look for block matches and mismatches in MD:Z string
            mdrun = 0
            runleft = run
            while runleft > 0 and mdo < len(mdp):
                op_m, run_m, st_m = mdp[mdo]
                run_comb = min(runleft, run_m)
                runleft -= run_comb
                assert op_m == 0 or op_m == 1
                rds += seq[rdoff:rdoff + run_comb]
                if op_m == 0:
                    rfs += seq[rdoff:rdoff + run_comb]
                else:
                    assert len(st_m) == run_comb
                    rfs += st_m
                    for i in xrange(0, run_comb):
                        assert st_m[i] != seq[rdoff+i:rdoff+i+1]
                mdrun += run_comb
                rdoff += run_comb
                # A stretch of matches in MD:Z could span M and I sections of
                # CIGAR
                if run_comb < run_m:
                    assert op_m == 0
                    mdp[mdo][1] -= run_comb
                else:
                    mdo += 1
        elif op == 1: # I
            rds += seq[rdoff:rdoff + run]
            rfs += "-" * run
            rdoff += run
        elif op == 2: # D
            op_m, run_m, st_m = mdp[mdo]
            assert op_m == 2
            assert run == run_m
            assert len(st_m) == run
            mdo += 1
            rds += "-" * run
            rfs += st_m
        elif op == 3: # N
            # If this is really long, this probably isn't the best way to do
            # this
            rds += "-" * run
            rfs += "-" * run
        elif op == 4: # S
            rdoff += run
        elif op == 5: # H
            pass
        elif op == 6: # P
            assert False
        elif op == 7: # =
            assert False
        elif op == 8: # X
            assert False
        else: assert False
    assert mdo == len(mdp)
    return rds, rfs

class ReservoirSampler(object):
    """ Simple reservoir sampler """
    def __init__(self, k):
        self.k = k
        self.r = []
        self.n = 0
    
    def add(self, obj):
        if self.n < self.k:
            self.r.append(obj)
        else:
            j = random.randint(0, self.n)
            if j < self.k:
                self.r[j] = obj
        self.n += 1
    
    def draw(self):
        return random.choice(self.r)
    
    def __len__(self):
        return self.n

class ScoreDist(object):
    """ Capture a map from scores observed to patterns of mismatches
        and gaps observed.  We do this instead of trying to generate
        patterns of gaps and mismatches to achieve a desired score,
        because the latter has many issues.  Chief amongst them:
        1. Hard to  """
    
    def __init__(self, k=10000):
        self.res = ReservoirSampler(k)
    
    def draw(self):
        assert len(self.res) > 0
        return self.res.draw()
    
    def add(self, al):
        # Extract quality string
        assert al.cigar is not None
        assert al.mdz is not None
        # Extract CIGAR, MD:Z
        cigarList = cigarToList(al.cigar)
        mdzList = mdzToList(al.mdz)
        sc = al.bestScore
        # Get stacked alignment
        rdAln, rfAln = cigarMdzToStacked(al.seq, cigarList, mdzList)
        self.res.add((al.fw, al.qual, rdAln, rfAln, sc))

class ScoreDistPaired(object):
    """ Capture a map from scores observed to patterns of mismatches
        and gaps observed.  We do this instead of trying to generate
        patterns of gaps and mismatches to achieve a desired score,
        because the latter has many issues.  Chief amongst them:
        1. Hard to  """
    
    def __init__(self, k=10000):
        self.res = ReservoirSampler(k)
    
    def draw(self):
        assert len(self.res) > 0
        return self.res.draw()
    
    def add(self, al1, al2):
        # Extract quality string
        assert al1.cigar is not None and al2.cigar is not None
        assert al1.mdz is not None and al2.mdz is not None
        def aln(al):
            cigarList = cigarToList(al.cigar)
            mdzList = mdzToList(al.mdz)
            return cigarMdzToStacked(al.seq, cigarList, mdzList)
        rdAln1, rfAln1 = aln(al1)
        rdAln2, rfAln2 = aln(al2)
        self.res.add((abs(al1.tlen),
                      (al1.fw, al1.qual, rdAln1, rfAln1, al1.bestScore),
                      (al2.fw, al2.qual, rdAln2, rfAln2, al2.bestScore)))

class Dist(object):
    """ Capture an empirical distribution.  Basically a histogram. """
    
    def __init__(self):
        self.hist = dict()
        self.tot = 0
        self.changed = True
        self.gen = None
        self.vals = None
    
    def __str__(self):
        return str(self.hist)
    
    def __len__(self):
        return self.tot
    
    def draw(self):
        if len(self) == 0:
            raise RuntimeError("Attempt to draw from empty empirical distribution")
        if self.changed:
            self.gen = WeightedRandomGenerator(self.hist.itervalues())
            self.changed = False
        return self.hist.keys()[self.gen.next()]
    
    def add(self, key):
        self.hist[key] = self.hist.get(key, 0) + 1
        self.tot += 1
        self.changed = True

class SequenceSimulator(object):
    """ Class that, given a collection of FASTA files, samples intervals of
        specified length from the strings contained in them. """
    
    def __init__(self, fafns, pickleFn, idx_fn=None, verbose=False):
        self.__re = re.compile('[^ACGTacgt]')
        self.refs = dict()
        self.names = []
        self.lens = []
        pickleExists = False
        if pickleFn is not None:
            pickleExists = os.path.exists(pickleFn)
        if pickleFn is not None and pickleExists:
            self.refs, self.names, self.lens = cPickle.load(open(pickleFn, 'rb'))
        else:
            totlen = 0
            pt_sz = 50000000
            last_pt = 0
            max_bases = sys.maxint
            if args.max_bases is not None:
                max_bases = args.max_bases
            abort = False
            for fafn in fafns:
                fafh = open(fafn, 'r')
                name = None
                for line in fafh:
                    line = line.rstrip()
                    ln = len(line)
                    if ln > 0 and line[0] == '>':
                        ind = line.find(" ")
                        if ind == -1: ind = len(line)
                        line = line[1:ind]
                        name = line
                        self.refs[name] = []
                        self.names.append(name)
                        self.lens.append(0)
                    else:
                        assert name is not None
                        self.refs[name].append(line)
                        self.lens[-1] += ln
                        totlen += ln
                        if verbose:
                            pt = totlen / pt_sz
                            if pt > last_pt:
                                print >> sys.stderr, "Read %d FASTA bytes..." % totlen
                            last_pt = pt
                        if totlen > max_bases:
                            abort = True
                            break
                fafh.close()
                if abort:
                    break
            for k in self.refs.iterkeys():
                self.refs[k] = ''.join(self.refs[k])
        if pickleFn is not None and not pickleExists:
            cPickle.dump((self.refs, self.names, self.lens), open(pickleFn, 'wb'), cPickle.HIGHEST_PROTOCOL)
        self.rnd = WeightedRandomGenerator(self.lens)
    
    def sim(self, ln, verbose=False):
        if verbose: print >>sys.stderr, "sim called..."
        # Pick a reference sequence in a weighted random fashon
        refi = self.rnd.next()
        assert refi < len(self.names)
        fw = True
        nm = self.names[refi] # reference name
        refoff = random.randint(0, self.lens[refi] - ln) # pick offset
        seq = self.refs[nm][refoff:refoff+ln] # extract substring
        # Simulated read can't overlap non-A-C-G-T character in reference
        while self.__re.search(seq):
            refoff = random.randint(0, self.lens[refi] - ln) # pick new offset
            seq = self.refs[nm][refoff:refoff+ln] # extract substring again
        seq = seq.upper()
        if random.random() > 0.5: # possibly reverse-complement
            fw = False
            seq = revcomp(seq) # reverse complement
        assert not "N" in seq
        if verbose:
            print >>sys.stderr, "...done"
        return (nm, refoff, fw, seq) # return ref id, ref offset, orientation, sequence

def mutate(rd, rdfw, scDistDraw):
    """ Given a read that already has the appropriate length (i.e.
        equal to the number of characters on the reference side of
        the alignment) 
        
        Modifies 'rd' in place. """ 
    fw, qual, rdAln, rfAln, sc = scDistDraw
    assert "N" not in rd.seq
    assert len(rd.seq) == len(rfAln) - rfAln.count('-')
    if rdfw != fw:
        qual, rdAln, rfAln = qual[::-1], rdAln[::-1], rfAln[::-1]
    rd.qual = qual # Install qualities
    # Walk along stacked alignment, making corresponding changes to
    # read sequence
    seq = []
    i, rdi, rfi = 0, 0, 0
    while i < len(rdAln):
        if rdAln[i] == '-':
            rfi += 1
        elif rfAln[i] == '-':
            seq.append(rdAln[i])
            rdi += 1
        elif rdAln[i] != rfAln[i] and rdAln[i] == 'N':
            seq.append('N')
            rfi += 1; rdi += 1
        elif rdAln[i] != rfAln[i]:
            assert rfi < len(rd.seq)
            oldc = rd.seq[rfi].upper()
            cs = ['A', 'C', 'G', 'T']
            assert oldc in cs, "oldc was %s" % oldc
            cs.remove(oldc)
            newc = random.choice(cs)
            seq.append(newc)
            rfi += 1; rdi += 1
        else:
            assert rfi < len(rd.seq)
            seq.append(rd.seq[rfi])
            rfi += 1; rdi += 1
        i += 1
    rd.seq = ''.join(seq)
    assert len(rd.seq) == len(rd.qual)

class SimulatorWrapper(object):
    
    """ Wrapper that sends requests to the Simualtor but uses information
        gathered during alignment so far to select such parameters as read
        length, concordant/discordant fragment length, etc. """
    
    def __init__(self, sim, m1fw, m2fw, tyd, scd, scd1, scd2, scdp, fld):
        self.sim  = sim  # sequence simulator
        self.m1fw = m1fw # whether mate 1 is revcomped w/r/t fragment
        self.m2fw = m2fw # whether mate 2 is revcomped w/r/t fragment
        self.tyd  = tyd  # type distribution (UU/CP/DP/UP)
        self.scd  = scd  # qualities/edits for unpaired reads
        self.scd1 = scd1 # qualities/edits for discordant mate 1
        self.scd2 = scd2 # qualities/edits for discordant mate 2
        self.scdp = scdp # qualities/edits for concordant pairs
        self.fld  = fld  # fragment length distribution
    
    def next(self):
        """ Simulate the next read/pair and associated quality values.  Return
            the simulated read along with information about where it
            originated. """
        ty = self.tyd.draw()
        if ty[1] == 'U':
            # Simulating unpaired read
            scDraw = self.scd.draw()
            _, _, _, rfAln, _ = scDraw
            rl = len(rfAln) - rfAln.count('-')
            refid, refoff, fw, seq = self.sim.sim(rl) # simulate it
            assert rl == len(seq)
            _, _, _, _, sc = scDraw
            read = Read.fromSimulator(seq, None, refid, refoff, fw, sc)
            mutate(read, fw, scDraw) # mutate unpaired read
            assert read.qual is not None
            return read, None
        else:
            # Simulating paired-end read
            fl, sc1Draw, sc2Draw = self.scdp.draw()
            _, _, _, rfAln1, _ = sc1Draw
            _, _, _, rfAln2, _ = sc2Draw
            rl1 = len(rfAln1) - rfAln1.count('-')
            rl2 = len(rfAln2) - rfAln2.count('-')
            refid, refoff, fw, seq = self.sim.sim(fl) # simulate fragment
            assert len(seq) == fl
            assert fl >= rl1
            assert fl >= rl2
            # get mates from fragment
            seq1, seq2 = seq[:rl1], seq[-rl2:]
            if not self.m1fw: seq1 = revcomp(seq1)
            if not self.m2fw: seq2 = revcomp(seq2)
            # Now we have the Watson offset for one mate or the other,
            # depending on which mate is to the left w/r/t Watson.
            refoff1, refoff2 = refoff, refoff
            if fw: refoff2 = refoff + fl - rl2
            else:  refoff1 = refoff + fl - rl1
            _, _, _, _, sc1 = sc1Draw
            _, _, _, _, sc2 = sc2Draw
            rd1 = Read.fromSimulator(seq1, None, refid, refoff1, fw, sc1)
            rd2 = Read.fromSimulator(seq2, None, refid, refoff2, fw, sc2)
            mutate(rd1, fw == self.m1fw, sc1Draw) # mutate mate 1
            mutate(rd2, fw == self.m2fw, sc2Draw) # mutate mate 2
            assert rd1.qual is not None
            assert rd2.qual is not None
            return rd1, rd2

class Input(object):
    """ Class that parses reads from input files and yields the reads/pairs
        produced using generators """
    
    @staticmethod
    def fastaParse(fh):
        lns = [ fh.readline().rstrip() for x in xrange(0, 2) ]
        orig = '\n'.join(lns) + '\n'
        if len(lns[0]) == 0: return None
        return Read(lns[0][1:], lns[1], None, orig)

    @staticmethod
    def fastqParse(fh):
        lns = [ fh.readline().rstrip() for x in xrange(0, 4) ]
        orig = '\n'.join(lns) + '\n'
        if len(lns[0]) == 0: return None
        return Read(lns[0][1:], lns[1], lns[3], orig)
    
    def __init__(self, format="fastq", unpFns=None, m1Fns=None, m2Fns=None):
        self.format = format
        if format == "fastq":
            self.parse = self.fastqParse
        elif format == "fasta":
            self.parse = self.fastaParse
        else:
            raise RuntimeError("Bad input format: '%s'" % format)
        self.unpFns = unpFns
        self.m1Fns = m1Fns
        self.m2Fns = m2Fns
    
    def __iter__(self):
        # Yield all the unpaired reads first
        if self.unpFns is not None:
            for unpFn in self.unpFns:
                with open(unpFn, 'r') as unpFh:
                    while True:
                        rd = self.parse(unpFh)
                        if rd is not None: yield (rd, None)
                        else: break # next file
        # Yield all the paired-end reads
        if self.m1Fns is not None:
            assert self.m2Fns is not None
            for (m1Fn, m2Fn) in zip(self.m1Fns, self.m2Fns):
                with open(m1Fn, 'r') as m1Fh:
                    with open(m2Fn, 'r') as m2Fh:
                        while True:
                            rd1, rd2 = self.parse(m1Fh), self.parse(m2Fh)
                            if rd1 is not None: yield (rd1, rd2)
                            else: break # next pair of files

class InputWrapper(object):
    """ Wraps the input reader so that we can eavesdrop on the reads and use
        them to build empirical distributions.  Also allows us to inject
        simulated reads now and again. """
    def __init__(self, inp, rddistUnp, rddistM1, rddistM2, qdistUnp, qdistM1, qdistM2):
        self.inp = inp
        self.rddistUnp, self.rddistM1, self.rddistM2 = rddistUnp, rddistM1, rddistM2
        self.qdistUnp, self.qdistM1, self.qdistM2 = qdistUnp, qdistM1, qdistM2
    
    def __iter__(self):
        for (rd1, rd2) in self.inp:
            assert rd1 is not None
            # ignore simulated reads
            if not rd1.name.startswith('!!ts!!'):
                if rd2 is not None:
                    # Paired-end case
                    # add read lengths to empirical distributions
                    self.rddistM1.add(len(rd1))
                    self.rddistM2.add(len(rd2))
                    # add quality values to empirical distributions
                    for q in map(lambda x: ord(x)-33, rd1.qual): self.qdistM1.add(q)
                    for q in map(lambda x: ord(x)-33, rd2.qual): self.qdistM2.add(q)
                else:
                    # Unpaired case
                    # add read length to empirical distribution
                    self.rddistUnp.add(len(rd1))
                    # add quality values to empirical distribution
                    for q in map(lambda x: ord(x)-33, rd1.qual): self.qdistUnp.add(q)
            yield (rd1, rd2)

class KNN(object):
    """ Encapsulates a K-nearest-neighbor classifier for estimating
        probability that an alignment tuple corresponds to a correct
        alignment.  Maintains an LRU cache mapping tuples to responses
        so we don't necessarily have to do the nearest-neighbor query
        each time.
        
        Caching code borrowed from:
        http://code.activestate.com/recipes/577970-simplified-lru-cache/
    """ 
    
    def __init__(self, n_neighbors, weights, maxsize=2048, sanity=False):
        self.classifier = KNeighborsClassifier(\
            n_neighbors=n_neighbors, warn_on_equidistant=False, weights=weights)
        
        self.didFit = False
        self.hits = 0
        self.misses = 0
        self.prev, self.prevProb = None, None
        self.sanity = sanity # do sanity-checking on cache results
        
        # Cache is a double-linked list
        # Link layout:     [PREV, NEXT, KEY, RESULT]
        self.root = root = [None, None, None, None]
        self.cache = cache = {}
        last = root
        for i in range(maxsize):
            key = object()
            cache[key] = last[1] = last = [last, root, key, None]
        root[0] = last
    
    def probCorrect(self, recl):
        """ Return the probability that the test tuple 'rec'
            corresponds to a correct alignment """
        assert self.didFit
        if self.prev is not None and recl == self.prev:
            return self.prevProb
        self.prev = recl
        rect = tuple(recl)
        
        cache = self.cache
        root = self.root
        link = cache.get(rect)
        if link is not None:
            # Cache hit
            link_prev, link_next, _, result = link
            link_prev[1] = link_next
            link_next[0] = link_prev
            last = root[0]
            last[1] = root[0] = link
            link[0] = last
            link[1] = root
            self.hits += 1
            self.prevProb = result
            if self.sanity:
                assert result == self.classifier.predict_proba([recl])[0][-1]
            return result
        # Cache miss
        result = self.classifier.predict_proba([recl])[0][-1]
        root[2] = rect
        root[3] = result
        oldroot = root
        root = self.root = root[1]
        root[2], oldkey = None, root[2]
        root[3], oldvalue = None, root[3]
        del cache[oldkey]
        cache[rect] = oldroot
        self.misses += 1
        self.prevProb = result
        return result
    
    def fit(self, train, lab):
        """ Fit the model (build the KDTree) """
        self.didFit = True
        self.classifier.fit(train, lab)

class Training(object):
    
    """ Encapsulates all the training data and all the classifiers. """
    def __init__(self):
        # These tables are for collecting training data pertaining to
        # individual reads and mates
        self.trainUnp,    self.labUnp,    self.classUnp    = [], [], None
        self.trainM1Disc, self.labM1Disc, self.classM1Disc = [], [], None
        self.trainM2Disc, self.labM2Disc, self.classM2Disc = [], [], None
        self.trainM1Conc, self.labM1Conc, self.classM1Conc = [], [], None
        self.trainM2Conc, self.labM2Conc, self.classM2Conc = [], [], None
        # Following tables are for a second layer of training for
        # concordantly-aligned pairs
        self.trainConc,   self.labConcM1,   self.labConcM2 = [], [], []
        self.classConcM1, self.classConcM2 = None, None
        self.trainConcFraglen = []
        # These scale factors are set for real when we fit
        self.scaleAs, self.scaleDiff = 1.0, 1.0
    
    def __len__(self):
        """ Return number of pieces of training data added so far """
        return len(self.trainUnp) + len(self.trainM1Disc) + len(self.trainM2Disc) + \
            len(self.trainM1Conc) + len(self.trainM2Conc)
    
    def add(self, al, correct):
        """ Add an alignment for a simulated read to our collection of
            training data. """
        rec = TrainingRecord.fromAlignment(al)
        if al.concordant:
            if al.mate1:
                self.trainM1Conc.append(rec.toList())
                self.labM1Conc.append(correct)
                self.trainConcFraglen.append(abs(al.tlen))
            else:
                self.trainM2Conc.append(rec.toList())
                self.labM2Conc.append(correct)
        elif al.discordant or al.unpPair:
            if al.mate1:
                self.trainM1Disc.append(rec.toList())
                self.labM1Disc.append(correct)
            else:
                self.trainM2Disc.append(rec.toList())
                self.labM2Disc.append(correct)
        else:
            self.trainUnp.append(rec.toList())
            self.labUnp.append(correct)
    
    def probCorrect(self, al1, al2=None):
        """ Return probability that given alignment is """
        assert al1.isAligned()
        rec = TrainingRecord.fromAlignment(al1, self.scaleAs, self.scaleDiff)
        if al2 is not None:
            assert al2.isAligned()
            assert al1.concordant
            assert al2.concordant
            assert al1.mate1 and not al1.mate2
            assert al2.mate2 and not al2.mate1
            rec1, rec2 = rec, TrainingRecord.fromAlignment(al2, self.scaleAs, self.scaleDiff)
            probCor1 = self.classM1Conc.probCorrect(rec1.toList())
            probCor2 = self.classM2Conc.probCorrect(rec2.toList())
            concRec = [ probCor1, probCor2, abs(al1.tlen) ]
            newProbCor1 = self.classConcM1.probCorrect(concRec)
            newProbCor2 = self.classConcM2.probCorrect(concRec)
            return newProbCor1, newProbCor2
        elif al1.discordant or al1.unpPair:
            classifier = self.classM1Disc if al1.mate1 else self.classM2Disc
            return classifier.probCorrect(rec.toList())
        else:
            return self.classUnp.probCorrect(rec.toList())
    
    def fit(self, num_neighbors, weights, scaleAs=1.0, scaleDiff=1.0):
        """ Train our KNN classifiers """
        self.scaleAs = scaleAs
        self.scaleDiff = scaleDiff
        
        def __adjustScale(t):
            # Scale fields as requested
            if scaleAs != 1.0 or scaleDiff != 1.0:
                for i in xrange(0, len(t)):
                    t[i][1] *= scaleAs
                    t[i][2] *= scaleDiff
        
        # Create and train each classifier
        if len(self.trainUnp) > 0:
            self.classUnp = KNN(n_neighbors=num_neighbors, weights=weights)
            __adjustScale(self.trainUnp)
            self.classUnp.fit(self.trainUnp, self.labUnp)
        if len(self.trainM1Disc) > 0:
            self.classM1Disc = KNN(n_neighbors=num_neighbors, weights=weights)
            __adjustScale(self.trainM1Disc)
            self.classM1Disc.fit(self.trainM1Disc, self.labM1Disc)
        if len(self.trainM2Disc) > 0:
            self.classM2Disc = KNN(n_neighbors=num_neighbors, weights=weights)
            __adjustScale(self.trainM2Disc)
            self.classM2Disc.fit(self.trainM2Disc, self.labM2Disc)
        if len(self.trainM1Conc) > 0:
            self.classM1Conc = KNN(n_neighbors=num_neighbors, weights=weights)
            __adjustScale(self.trainM1Conc)
            self.classM1Conc.fit(self.trainM1Conc, self.labM1Conc)
        if len(self.trainM2Conc) > 0:
            self.classM2Conc = KNN(n_neighbors=num_neighbors, weights=weights)
            __adjustScale(self.trainM2Conc)
            self.classM2Conc.fit(self.trainM2Conc, self.labM2Conc)
        
        assert len(self.trainM1Conc) == len(self.trainM2Conc)
        assert len(self.trainM1Conc) == len(self.trainConcFraglen)
        assert len(self.trainM1Disc) == len(self.trainM2Disc)
        
        # Create training data for the concordant-alignment classifier.
        # This depends on M1Conc and M2Conc classifiers already having
        # been fit.
        if len(self.trainM1Conc) > 0:
            for i in xrange(0, len(self.trainM1Conc)):
                m1c, m2c = self.trainM1Conc[i], self.trainM2Conc[i]
                p1 = self.classM1Conc.probCorrect(m1c)
                p2 = self.classM2Conc.probCorrect(m2c)
                self.trainConc.append((p1, p2, self.trainConcFraglen[i]))
                self.labConcM1.append(self.labM1Conc[i])
                self.labConcM2.append(self.labM2Conc[i])
        
            self.classConcM1 = KNN(n_neighbors=num_neighbors, weights=weights)
            self.classConcM2 = KNN(n_neighbors=num_neighbors, weights=weights)
            
            self.classConcM1.fit(self.trainConc, self.labM1Conc)
            self.classConcM2.fit(self.trainConc, self.labM2Conc)
    
    def save(self, fn):
        """ Save all training data to a file """
        save = (\
        self.trainUnp,    self.labUnp, \
        self.trainM1Disc, self.labM1Disc, \
        self.trainM2Disc, self.labM2Disc, \
        self.trainM1Conc, self.labM1Conc, \
        self.trainM2Conc, self.labM2Conc, \
        self.trainConc,   self.labConcM1,   self.labConcM2, \
        self.classConcM1, self.classConcM2, \
        self.trainConcFraglen)
        with open(fn, 'wb') as trainOfh: cPickle.dump(save, trainOfh, cPickle.HIGHEST_PROTOCOL)
    
    def load(self, fn):
        """ Load all training data from a file """
        with open(fn, 'rb') as trainIfh:
            (\
            self.trainUnp,    self.labUnp, \
            self.trainM1Disc, self.labM1Disc, \
            self.trainM2Disc, self.labM2Disc, \
            self.trainM1Conc, self.labM1Conc, \
            self.trainM2Conc, self.labM2Conc, \
            self.trainConc,   self.labConcM1,   self.labConcM2, \
            self.classConcM1, self.classConcM2, \
            self.trainConcFraglen) = cPickle.load(trainIfh)
    
    def hits(self):
        return self.classConcM1.hits if self.classConcM1 is not None else 0 + \
               self.classConcM2.hits if self.classConcM2 is not None else 0 + \
               self.classUnp.hits    if self.classUnp    is not None else 0 + \
               self.classM1Disc.hits if self.classM1Disc is not None else 0 + \
               self.classM2Disc.hits if self.classM2Disc is not None else 0 + \
               self.classM1Conc.hits if self.classM1Conc is not None else 0 + \
               self.classM2Conc.hits if self.classM2Disc is not None else 0
    
    def misses(self):
        return self.classConcM1.misses if self.classConcM1 is not None else 0 + \
               self.classConcM2.misses if self.classConcM2 is not None else 0 + \
               self.classUnp.misses    if self.classUnp    is not None else 0 + \
               self.classM1Disc.misses if self.classM1Disc is not None else 0 + \
               self.classM2Disc.misses if self.classM2Disc is not None else 0 + \
               self.classM1Conc.misses if self.classM1Conc is not None else 0 + \
               self.classM2Conc.misses if self.classM2Disc is not None else 0

class TrainingRecord(object):
    """ A single tuple of per-read training data including:
        1. Read length
        . Read's best score
        34. Read's second-best score """
    def __init__(self, rdlen, bestSc, secbestSc, scaleAs=1.0, scaleDiff=1.0):
        self.rdlen = rdlen
        self.bestSc = float(bestSc) * scaleAs
        self.secbestSc = secbestSc
        if secbestSc is not None:
            self.scDiff = float(abs(self.bestSc - self.secbestSc))
        else:
            self.scDiff = 10000.0
        self.scDiff *= scaleDiff
        assert self.repOk()
    
    @classmethod
    def fromAlignment(cls, al, scaleAs=1.0, scaleDiff=1.0):
        """ Initialize training record with respect to an alignment and a
            boolean indicating whether it is correct. """
        secbest = al.secondBestScore
        if secbest is None or al.thirdBestScore > secbest:
            secbest = al.thirdBestScore
        return cls(len(al), al.bestScore, secbest,
                   scaleAs=scaleAs, scaleDiff=scaleDiff)
    
    def toList(self):
        """ Return simple list form """
        return [ self.rdlen, self.bestSc, self.scDiff ]
    
    def repOk(self):
        """ Check for internal consistency """
        assert self.rdlen is not None
        assert self.bestSc is not None
        assert self.scDiff is not None
        return True

class Output(threading.Thread):
    """ Encapsulates the output reader.  Reads SAM output from the aligner,
        updates empirical distributions for e.g. fragment length, and looks for
        records that correspond to simulated reads.  If a record corresponds to
        a simulated read, its correctness will be checked and a tuple
        written """
    
    def __init__(self, samIfh, samOfh, unalFh, trainSink,
                 scDistUnp, scDistM1, scDistM2, scDistPair, fragDist, typeDist):
        threading.Thread.__init__(self)
        self.samIfh = samIfh         # SAM records come from here
        self.samOfh = samOfh         # normal (non-simulated) SAM recs go here
        self.trainSink = trainSink   # write training data here
        self.fragDist = fragDist     # fragment length distribution
        self.scDistUnp = scDistUnp   # qualities/edits for unpaired
        self.scDistM1 = scDistM1     # qualities/edits for mate #1s
        self.scDistM2 = scDistM2     # qualities/edits for mate #2s
        self.scDistPair = scDistPair # qualities/edits for concordant pairs
        self.typeDist = typeDist     # alignment type (UU/CP/DP/UP)
        self.unalFh = unalFh         # write unaligned reads here
        self.scDiffs = {}
    
    def run(self):
        lastAl = None
        for ln in self.samIfh:
            if ln[0] == '@':
                if self.samOfh is not None:
                    self.samOfh.write(ln) # header line
                continue
            al = Alignment(ln)
            nm, flags, refid, pos, _, _, _, _, _, seq, qual, _ = string.split(ln, '\t', 11)
            if al.name[0] == '!' and al.name.startswith('!!ts!!'):
                # this is a simulated read
                _, refid, fw, refoff, sc = string.split(al.name, '!!ts-sep!!')
                sc = int(sc)
                refoff = int(refoff)
                if al.isAligned():
                    scDiff = sc - al.bestScore
                    self.scDiffs[scDiff] = self.scDiffs.get(scDiff, 0) + 1
                    correct = False
                    # Check reference id, orientation
                    if refid == al.refid and fw == al.orientation():
                        # Check offset
                        correct = abs(refoff - al.pos) < args.wiggle
                    if self.trainSink is not None:
                        self.trainSink(al, correct)
                elif self.unalFh is not None:
                    # Perhaps write unaligned simulated read to file
                    self.unalFh.write("@%s\n%s\n+\n%s\n" % (nm, seq, qual))
            else:
                # Take alignment info into account
                if al.isAligned():
                    if lastAl is not None and al.name == lastAl.name and al.alType == "CP":
                        assert lastAl.alType == "CP"
                        mate1, mate2 = al, lastAl
                        if (lastAl.flags & 64) != 0:
                            mate1, mate2 = mate2, mate1
                        self.fragDist.add(abs(mate1.tlen)) # alignment pair
                        self.scDistPair.add(mate1, mate2)
                    elif (al.flags & 64) != 0:
                        self.scDistM1.add(al)
                    elif (al.flags & 128) != 0:
                        self.scDistM2.add(al)
                    else:
                        self.scDistUnp.add(al)
                    self.typeDist.add(al.alType)
                # Send SAM to SAM output filehandle
                if self.samOfh is not None:
                    self.samOfh.write(ln)
            lastAl = al

class AsyncWriter(threading.Thread):
    def __init__(self, fh, q, name=""):
        threading.Thread.__init__(self)
        self.fh = fh
        self.q = q
        self.name = name
    
    def run(self):
        i = 0
        while True:
            item = self.q.get()
            if item is None:
                self.q.task_done()
                self.fh.flush()
                self.fh.close()
                break
            self.fh.write(item)
            self.fh.flush()
            self.q.task_done()
            i += 1

def createInput(rddistUnp, rddistM1, rddistM2, qdistUnp, qdistM1, qdistM2):
    """ Return an Input object that reads all user-provided input reads """
    if args.fastq:
        return InputWrapper(\
            Input(format="fastq", unpFns=args.U, m1Fns=args.m1, m2Fns=args.m2),
            rddistUnp, rddistM1, rddistM2, qdistUnp, qdistM1, qdistM2)
    elif args.fasta:
        return InputWrapper(\
            Input(format="fasta", unpFns=args.U, m1Fns=args.m1, m2Fns=args.m2),
            rddistUnp, rddistM1, rddistM2, qdistUnp, qdistM1, qdistM2)

def go():
    """ Main driver for tandem simulator """
    
    import tempfile
    
    if args.ref is None and args.training_input is None:
        raise RuntimeError("Must specify --ref")
    
    random.seed(args.seed)
    
    # When do we start generating training data?  If we train as we go,
    # we can start refining MAPQ estimates during the first pass, but we have
    # less data to direct the simulator to the most appropriate portions of the
    # training space.  If we wait until all the reads have been aligned and
    # generate it then, we need a second pass to refine any/all of the MAPQ
    # estimates but we have a lot of data to direct the simulator to the most
    # relevant portions of the training space. 
    
    # Construct command to invoke bowtie2
    bt2_cmd = "bowtie2 "
    if args.bt2_exe is not None:
        bt2_cmd = args.bt2_exe + " "
    if args.bt2_args is not None:
        bt2_cmd += args.bt2_args
    bt2_cmd += " --reorder --sam-no-qname-trunc -q "
    
    fifoArgs = ""
    fifoFns = [None, None, None] # FIFO filenames
    fifoFhs = [None, None, None] # FIFO filehandles
    fifoQs  = [None, None, None] # FIFO read queues
    fifoWs  = [None, None, None] # FIFO worker threads
    
    def makeFifos():
        fifoArgs = ''
        # Make temporary directory to store FIFOs
        tmpdir = tempfile.mkdtemp()
        # Make the appropriate FIFOs given the types of input we have
        if args.U is not None:
            fifoFns[0] = os.path.join(tmpdir, 'U.fifo')
            os.mkfifo(fifoFns[0])
            fifoArgs += (" -U " + fifoFns[0])
        if args.m1 is not None:
            fifoFns[1] = os.path.join(tmpdir, 'm1.fifo')
            os.mkfifo(fifoFns[1])
            fifoArgs += (" -1 " + fifoFns[1])
        if args.m2 is not None:
            fifoFns[2] = os.path.join(tmpdir, 'm2.fifo')
            os.mkfifo(fifoFns[2])
            fifoArgs += (" -2 " + fifoFns[2])
        return fifoArgs
    
    def openBowtie(bt2_cmd):
        bt2_cmd += makeFifos()
        # Open the Bowtie 2 process, which is going to want to start reading from
        # one or more of the FIFOs
        pipe = subprocess.Popen(bt2_cmd, shell=True, stdout=subprocess.PIPE)
        if args.U  is not None:
            fifoFhs[0] = open(fifoFns[0], 'w')
        if args.m1 is not None:
            fifoFhs[1] = open(fifoFns[1], 'w')
        if args.m2 is not None:
            fifoFhs[2] = open(fifoFns[2], 'w')
        # For each input type (unpaired, mate1, mate2), initialize a queue and a
        # thread that takes reads from the queue and passes each along to the
        # appropriate FIFO.  It's important to have a separate thread for each FIFO
        # or we get deadlocks.
        for i in xrange(0, 3):
            if fifoFhs[i] is not None:
                fifoQs[i] = Queue()
                fifoWs[i] = AsyncWriter(fifoFhs[i], fifoQs[i], "Thread %d" % i)
                fifoWs[i].start()
        return pipe
    
    def closeFifos():
        # Write None to all the FIFOs to inform them we're done giving them reads
        for i in xrange(0, 3):
            if fifoFhs[i] is not None:
                fifoQs[i].put(None)
        # Join all the FIFOs and input handling threads, close and delete the FIFOs 
        for i in xrange(0, 3):
            if fifoFhs[i] is not None:
                fifoQs[i].join()
                if fifoWs[i] is not None:
                    fifoWs[i].join()
                    fifoWs[i] = None
                if fifoFhs[i] is not None:
                    fifoFhs[i].close()
                    fifoFhs[i] = None
                os.unlink(fifoFns[i])
    
    # Structure encapsulating training data and KNN classifiers
    training = Training()
    
    # Timing measurements
    setupIval, al1Ival, al2Ival, fitIval, samIval = None, None, None, None, None
    
    if args.training_input is None:
        
        # The training data is not being given as input, so we have to
        # generate it with two Bowtie 2 runs.  The first run collects
        # some informations about data distributions.  The second run
        # draws from those distributions to simulate reads.  Each
        # simulated read that aligns is used as a training tuple.
        
        st = time.clock()
        
        # ##################################################
        # ALIGN REAL DATA
        # ##################################################
        
        pipe = openBowtie(bt2_cmd)
        samIfh = tempfile.TemporaryFile()
        
        if args.verbose:
            print >> sys.stderr, "Real-data Bowtie 2 command: '%s'" % bt2_cmd
        
        scDistUnp, scDistM1, scDistM2, scDistPair, typeDist, fragDist = \
            ScoreDist(), ScoreDist(), ScoreDist(), ScoreDistPaired(), Dist(), Dist()
        
        # Create the thread that eavesdrops on output from bowtie2
        othread = Output(
            pipe.stdout,     # SAM input filehandle
            samIfh,          # SAM output filehandle
            None,            # Write unaligned reads here
            None,            # Training record sink
            scDistUnp,       # qualities/edits for unpaired
            scDistM1,        # qualities/edits for mate1
            scDistM2,        # qualities/edits for mate2
            scDistPair,      # qualities/edits for concordant pairs
            fragDist,        # Fragment dist
            typeDist)        # Alignment type dist
        othread.start()
        rddistUnp, rddistM1, rddistM2 = Dist(), Dist(), Dist()
        qdistUnp, qdistM1, qdistM2 = Dist(), Dist(), Dist()
        
        # Construct sequence and quality simulators
        sim = SequenceSimulator(args.ref, args.pickle_ref, idx_fn=args.ref_idx)
        simw = SimulatorWrapper(\
            sim,        # sequence simulator
            not args.m1rc, # whether m1 is revcomped w/r/t fragment
            args.m2fw,  # whether m2 is revcomped w/r/t fragment
            typeDist,   # alignment type distribution
            scDistUnp,  # qualities/edits for unpaired
            scDistM1,   # qualities/edits for mate 1
            scDistM2,   # qualities/edits for mate 2
            scDistPair, # qualities/edits for concordant pairs
            fragDist)   # fragment-length distribution
        
        # Stop timing setup
        setupIval = time.clock() - st
        st = time.clock()
        
        print >> sys.stderr, "Initializing threads, queues and FIFOs"
        
        # Read through all the input read files and direct all reads to the
        # appropriate queue
        upto = args.upto or sys.maxint
        numReads = 0
        for (rd1, rd2) in iter(createInput(rddistUnp, rddistM1, rddistM2, qdistUnp, qdistM1, qdistM2)):
            numReads += 1
            if rd2 is not None:
                # Write to -1/-2 filehandles
                assert fifoFhs[1] is not None
                assert fifoFhs[2] is not None
                fifoQs[1].put(str(rd1))
                fifoQs[2].put(str(rd2))
            else:
                # Write to -U filehandle
                assert fifoFhs[0] is not None
                fifoQs[0].put(str(rd1))
            if numReads >= upto:
                break
        
        print >> sys.stderr, "Finished reading real reads"
    
        closeFifos()
        othread.join() # join the thread that monitors aligner output
        print >> sys.stderr, "Finished closing data FIFOs and joining output-monitoring thread"
        
        # Stop timing interval for alignment phase 1
        al1Ival = time.clock() - st
        st = time.clock()
        
        # ##################################################
        # ALIGN SIMULATED DATA
        # ##################################################
        #
        # Now we re-open Bowtie 2 using the same arguments.  This time we only give
        # it simulated reads.
        #
        
        print >> sys.stderr, "Opening new Bowtie 2"
        pipe = openBowtie(bt2_cmd + " --mapq-extra")
        
        # Function gets called with each new piece of training data
        def trainingSink(al, correct): training.add(al, correct)
        
        unalFh = None
        if args.un_sim is not None:
            unalFh = open(args.un_sim, 'w')
        
        # Create the thread that eavesdrops on output from bowtie2 with simulated
        # input
        othread = Output(\
            pipe.stdout,     # SAM input filehandle
            None,            # SAM output filehandle
            unalFh,          # Write unaligned reads here
            trainingSink,    # Training record sink
            None,            # qualities/edits for unpaired
            None,            # qualities/edits for mate 1
            None,            # qualities/edits for mate 2
            None,            # qualities/edits for concordant pairs
            None,            # Fragment dist
            None)            # Alignment type dist
        othread.start()
        
        # Simulate reads from empirical distributions
        for i in xrange(0, args.num_reads):
            if (i+1 % 1000) == 0: print >> sys.stderr, "Generating read %d" % i
            rd1, rd2 = simw.next()
            if rd2 is not None:
                # Paired-end simulated read
                fifoQs[1].put(str(rd1))
                fifoQs[2].put(str(rd2))
            else:
                # Unpaired simulated read
                fifoQs[0].put(str(rd1))
        
        print >> sys.stderr, "Finished simulating reads"
        closeFifos()
        print >> sys.stderr, "Closed FIFOs"
        othread.join() # join the thread that monitors aligner output
        if unalFh is not None:
            unalFh.close()
        print >> sys.stderr, "Finished closing simulation FIFOs and joining output-monitoring thread"
        
        # Stop timing interval for alignment phase 2
        al2Ival = time.clock() - st
        
        print >> sys.stderr, "Score difference (expected - actual) histogram:"
        for k, v in sorted(othread.scDiffs.iteritems()):
            print >>sys.stderr, "  %d: %d" % (k, v)
        
        if args.save_training is not None: training.save(args.save_training)
    else:
        st = time.clock()
        
        # Training data is being given to us
        if args.sam_input is None:
            raise RuntimeError("Must specify --sam-input along with --training-input")
        training.load(args.training_input)
        print >> sys.stderr, "Read %d training tuples from '%s'" % (len(training), args.training_input)
        samIfh = open(args.sam_input, 'r')
        
        # Stop timing interval for setup phase
        setupIval = time.clock() - st
    
    # Build KNN classifiers
    st = time.clock()
    weights = 'distance' if args.distance_weight else 'uniform'
    training.fit(args.num_neighbors, weights, scaleAs=args.as_scale, scaleDiff=args.diff_scale)
    print >> sys.stderr, "Finished fitting KNN classifiers on %d training tuples" % len(training)
    fitIval = time.clock() - st
    st = time.clock()
    
    # TODO: Need to pair up alignments so we can potentially use the
    # concordant-alignment model to calculate their MAPQs
    nrecs = 0
    mapqDiff = 0.0
    npair, nunp = 0, 0
    
    with open(args.S, 'w') as samOfh: # open output file
        
        def emitNewSam(samrec, al, probCorrect):
            probIncorrect = 1.0 - probCorrect # convert to probability incorrect
            mapq = args.max_mapq
            if probIncorrect > 0.0:
                mapq = min(-10.0 * math.log10(probIncorrect), args.max_mapq)
            samrec = samrec.rstrip()
            xqIdx = samrec.find("\tXQ:f:")
            if xqIdx != -1:
                samrec = samrec[:xqIdx]
            samOfh.write("\t".join([samrec, "XQ:f:" + str(mapq)]) + "\n")
            return mapq - float(al.mapq)
        
        samIfh.seek(0)                # rewind temporary SAM file
        concMap = {}
        for samrec in samIfh:         # read each record
            if samrec[0] == '@':      # pass headers straight through
                samOfh.write(samrec)
                continue
            al = Alignment(samrec)    # parse alignment
            nrecs += 1
            if al.isAligned():
                # Is this alignment part of a concordantly-aligned paired-end read?
                probCorrect = None
                if al.concordant:
                    if al.name in concMap:
                        al1, samrec1 = al, samrec
                        al2, samrec2 = concMap[al.name]
                        npair += 2
                        del concMap[al.name]
                        if al2.mate1:
                            al1, al2 = al2, al1
                            samrec1, samrec2 = samrec2, samrec1
                        probCorrect1, probCorrect2 = training.probCorrect(al1, al2)
                        mapqDiff += emitNewSam(samrec1.rstrip(), al1, probCorrect1)
                        mapqDiff += emitNewSam(samrec2.rstrip(), al2, probCorrect2)
                    else:
                        concMap[al.name] = (al, samrec)
                else:
                    nunp += 1
                    probCorrect = training.probCorrect(al)
                    mapqDiff += emitNewSam(samrec.rstrip(), al, probCorrect)
            else:
                samOfh.write(samrec)
        assert len(concMap) == 0
    
    samIval = time.clock() - st
    
    print >> sys.stderr, "Finished writing final SAM output (%d records) to '%s'" % (nrecs, args.S)
    print >> sys.stderr, "  concordant-pair SAM records: %d" % npair
    print >> sys.stderr, "  non-concordant-pair SAM records: %d" % nunp
    print >> sys.stderr, "  total mapping quality difference (new - old): %0.3f" % mapqDiff
    print >> sys.stderr, "  total KNN hits=%d, misses=%d" % (training.hits(), training.misses()) 
    
    if setupIval is not None:
        print >> sys.stderr, "Setup running time: %0.3f secs" % setupIval 
    if al1Ival is not None:
        print >> sys.stderr, "Alignment (real reads): %0.3f secs" % al1Ival 
    if al2Ival is not None:
        print >> sys.stderr, "Alignment (simulated reads): %0.3f secs" % al2Ival 
    if fitIval is not None:
        print >> sys.stderr, "Fitting: %0.3f secs" % fitIval
    if samIval is not None:
        print >> sys.stderr, "MAPQ calculation / SAM rewriting: %0.3f secs" % samIval
    
    # Print timing results
    
    # Close temporary SAM output file; it will be deleted immediately
    samIfh.close()

if __name__ == "__main__":
    
    import argparse
    
    parser = argparse.ArgumentParser(\
        description='Evaluate the sensitivity of an alignment tool w/r/t given'
                    'genome and set of alignment parameters.')
    
    parser.add_argument(\
        '--ref', metavar='path', type=str, nargs='+', help='FASTA file(s) containing reference genome sequences')
    parser.add_argument(\
        '--pickle-ref', metavar='path', type=str,
        help='Pickle FASTA input for speed, or use pickled version if it '
             'exists already.  Pickled version is stored at given path')
    parser.add_argument(\
        '--ref-idx', metavar='path', type=str, help='File containing index for FASTA file(s)')
    parser.add_argument(\
        '--U', metavar='path', type=str, nargs='+', help='Unpaired read files')
    parser.add_argument(\
        '--S', metavar='path', type=str, required=True, help='Write SAM output here')
    parser.add_argument(\
        '--m1', metavar='path', type=str, nargs='+', help='Mate 1 files')
    parser.add_argument(\
        '--m2', metavar='path', type=str, nargs='+', help='Mate 2 files')
    parser.add_argument(\
        '--m1rc', action='store_const', const=True, default=False,
        help='Set if mate 1 is reverse-complemented w/r/t the fragment')
    parser.add_argument(\
        '--m2fw', action='store_const', const=True, default=False,
        help='Set if mate 2 is oriented forward w/r/t the fragment')
    parser.add_argument(\
        '--fasta', action='store_const', const=True, default=False, help='Reads are FASTA')
    parser.add_argument(\
        '--fastq', action='store_const', const=True, default=True, help='Reads are FASTQ')
    parser.add_argument(\
        '--max-mapq', metavar='float', type=float, default=100.0,
        required=False, help='Maximum MAPQ possible for an alignment')
    parser.add_argument(\
        '--seed', metavar='int', type=int, default=99099,
        required=False, help='Integer to initialize pseudo-random generator')
    parser.add_argument(\
        '--num-reads', metavar='int', type=int, default=100000,
        required=False, help='Number of reads to simulate')
    parser.add_argument(\
        '--num-neighbors', metavar='int', type=int, default=500,
        required=False, help='Number of neighbors to use for k-nearest-neighbors')
    parser.add_argument(\
        '--upto', metavar='int', type=int, default=None,
        required=False, help='Stop after this many input reads')
    parser.add_argument(\
        '--wiggle', metavar='int', type=int, default=30,
        required=False, help='Wiggle room to allow in starting position when determining whether alignment is correct')
    parser.add_argument(\
        '--max-ref-bases', metavar='int', dest='max_bases', type=int, default=None,
        required=False, help='Stop reading in FASTA once we exceed this many reference nucleotides')
    parser.add_argument(\
        '--training-input', metavar='path', type=str,
        help='Training data to use to predict new MAPQs for a SAM file.  Use with --sam-input.')
    parser.add_argument(\
        '--save-training', metavar='path', type=str,
        help='Save training data to file.')
    parser.add_argument(\
        '--sam-input', metavar='path', type=str,
        help='Input SAM file to apply training data to.  Use with --training-input.')
    parser.add_argument(\
        '--un-sim', metavar='path', type=str,
        help='Write unaligned simulated reads to this file.')
    parser.add_argument(\
        '--bt2-exe', metavar='path', dest='bt2_exe', type=str,
        help='Path to Bowtie 2 exe')
    parser.add_argument(\
        '--bt2-args', metavar='args', dest='bt2_args', type=str,
        help='Arguments to pass to Bowtie 2 (besides input an output)')
    parser.add_argument(\
        '--distance-weight', action='store_const', const=True, default=False,
        help='Do distance weighting when doing KNN')
    parser.add_argument(\
        '--as-scale', metavar='float', type=float, default=1.0,
        help='Multiply AS:i scores by this before making them into training records')
    parser.add_argument(\
        '--diff-scale', metavar='float', type=float, default=3.0,
        help='Multiply AS:i - XS:i scores by this before making them into training records')
    parser.add_argument(\
        '--sanity', dest='sanity', action='store_const', const=True, default=False,
        help='Do various sanity checks')
    parser.add_argument(\
        '--test', dest='test', action='store_const', const=True, default=False,
        help='Do unit tests')
    parser.add_argument(\
        '--profile', action='store_const', const=True,
        default=False, help='Print profiling info')
    parser.add_argument(\
        '--verbose', dest='verbose', action='store_const', const=True,
        default=False, help='Be talkative')
    parser.add_argument(\
        '--version', dest='version', action='store_const', const=True,
        default=False, help='Print version and quit')
    
    args = parser.parse_args()

    if args.test:
        
        import unittest
        
        class Test(unittest.TestCase):
            def test_1(self):
                dst = Dist()
                dst.add("CP")
                dst.add("CP")
                dst.add("DP")
                dr = dst.draw()
                self.assertTrue(dr == "CP" or dr == "DP")

        unittest.main(argv=[sys.argv[0]])
    elif args.profile:
        import cProfile
        cProfile.run('go()')
    else:
        go()
