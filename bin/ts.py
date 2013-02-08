#!/usr/bin/env python

"""
ts.py

A "tandem simulator," which wraps an alignment tool as it runs, eavesdrops on
the input and output, and builds a model that can be used to improve the
quality values calculated for aligned reads.

Right now we have to load the reference FASTA separately, and we need the user
to tell us about the scoring params so we can simulate reads with known best
score.

Things we learn from reads
==========================

- Read length distribution
-

Things we learn from alignments
===============================

- Fragment length distribution
- Score distribution (or # mismatch / # gap distribution)

"""

__author__ = "Ben Langmead"
__email__ = "langmea@cs.jhu.edu"

import os
import sys
import re
import argparse
import threading
import string
import subprocess
import tempfile
import signal
import traceback
from Queue import Queue

def quit_handler(signum,frame):
    traceback.print_stack()
signal.signal(signal.SIGQUIT,quit_handler)

parser = argparse.ArgumentParser(\
    description='Evaluate the sensitivity of an alignment tool w/r/t given'
                'genome and set of alignment parameters.')

parser.add_argument(\
    '--ref', metavar='path', type=str, nargs='+', required=True, help='FASTA file(s) containing reference genome sequences')
parser.add_argument(\
    '--U', metavar='path', type=str, nargs='+', help='Unpaired read files')
parser.add_argument(\
    '--S', metavar='path', type=str, help='Write SAM output here')
parser.add_argument(\
    '--m1', metavar='path', type=str, nargs='+', help='Mate 1 files')
parser.add_argument(\
    '--m2', metavar='path', type=str, nargs='+', help='Mate 2 files')
parser.add_argument(\
    '--fasta', action='store_const', const=True, default=False, help='Reads are FASTA')
parser.add_argument(\
    '--fastq', action='store_const', const=True, default=True, help='Reads are FASTQ')
parser.add_argument(\
    '--scoring', metavar='str', type=str, required=False,
    default='1,2,6,1,5,3,5,3',
    help='MatchBonus,MismatchMinPen,MismatchMaxPen,NPen,ReadGapConst,ReadGapLinear,RefGapConst,RefGapLinear')
parser.add_argument(\
    '--num-reads', metavar='int', dest='num_reads', type=int, default=100,
    required=False, help='Number of reads to simulate')
parser.add_argument(\
    '--max-ref-bases', metavar='int', dest='max_bases', type=int, default=None,
    required=False, help='Stop reading in FASTA once we exceed this many reference nucleotides')
parser.add_argument(\
    '--bt2-exe', metavar='path', dest='bt2_exe', type=str,
    help='Path to Bowtie 2 exe')
parser.add_argument(\
    '--bt2-args', metavar='args', dest='bt2_args', type=str,
    help='Arguments to pass to Bowtie 2 (besides input an output)')
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

class Read(object):
    """ Encapsulates one read """
    def __init__(self, name, seq, qual, orig=None):
        self.name = name
        self.seq = seq
        self.qual = qual
        self.orig = orig
    
    def __len__(self):
        """ Return number of nucleotides in read """
        return len(self.seq)
    
    def __str__(self):
        """ Return string representation """
        if self.orig is not None:
            return self.orig # original string preferred
        elif self.qual is not None:
            return "@%s\n%s\n+\n%s\n" % (self.name, self.seq, self.qual)
        else:
            return ">%s\n%s\n" % (self.name, self.seq)

class ReadPair(object):
    """ Encapsulates a pair of reads with mate1/mate2 labels """
    def __init__(self, rd1, rd2):
        self.rd1, self.rd2 = rd1, rd2

class Alignment(object):
    """ Encapsulates an alignment record for a single aligned read """

    __asRe = re.compile('AS:i:([-]?[0-9]+)')
    __xsRe = re.compile('XS:i:([-]?[0-9]+)')
    __ytRe = re.compile('YT:Z:([A-Z]+)')
    
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
        self.tlen = int(self.mapq)
        se = Alignment.__asRe.search(self.extra)
        self.bestScore = None
        if se is not None:
            self.bestScore = se.group(1)
        se = Alignment.__xsRe.search(self.extra)
        self.secondBestScore = None
        if se is not None:
            self.secondBestScore = se.group(1)
        se = Alignment.__ytRe.search(self.extra)
        self.alType = None
        if se is not None:
            self.alType = se.group(1)
        assert self.alType is not None
        assert not self.isAligned() or self.bestScore is not None, ln
    
    def isAligned(self):
        return (self.flags & 4) == 0
    
    def orientation(self):
        if (self.flags & 16) != 0:
            return "-"
        else:
            return "+"
    
    def __len__(self):
        return len(self.seq)

class AlignmentPair(object):
    """ Encapsulates a pair of alignments for two ends of a paired-end read """
    def __init__(self, al1, al2):
        self.al1, self.al2 = al1, al2
    
    def fraglen(self):
        return abs(self.al1.tlen)

class FragmentLengthDist(object):
    """ Capture empirical distribution of fragment lengths """
    
    def __init__(self):
        self.hist = dict()
    
    def __str__(self):
        return str(self.hist)
    
    def addFragment(self, alPair):
        """ Given an AlignmentPair, add the implied fragment length to our
            histogram """
        ln = alPair.fraglen()
        self.hist[ln] = self.hist.get(ln, 0) + 1

class ScoreDist(object):
    """ Capture empirical distribution of read lengths """
    
    def __init__(self):
        self.hist = dict()
    
    def __str__(self):
        return str(self.hist)
    
    def addAlignment(self, al):
        assert al.bestScore is not None
        sc = al.bestScore
        k = (sc, len(al))
        self.hist[k] = self.hist.get(k, 0) + 1

class ReadLengthDist(object):
    """ Capture empirical distribution of read lengths """
    
    def __init__(self):
        self.hist = dict()
    
    def __str__(self):
        return str(self.hist)
    
    def addRead(self, r):
        ln = len(r)
        self.hist[ln] = self.hist.get(ln, 0) + 1

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

class Simulator(object):
    
    """ Class that, given a collection of FASTA files, samples intervals of
        specified length from the strings contained in them. """
    
    def __init__(self, fafns, verbose=False):
        self.refs = dict()
        self.names = []
        self.lens = []
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
                if line.startswith('>'):
                    ind = line.index(" ")
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
        self.rnd = WeightedRandomGenerator(self.lens)
        self.__re = re.compile('[^ACGTacgt]')
    
    def sim(self, ln, verbose=False):
        if verbose:
            print >>sys.stderr, "sim called..."
        refi = self.rnd.next()
        assert refi < len(self.names)
        fw = True
        nm = self.names[refi]
        off = random.randint(0, self.lens[refi] - ln)
        seq = self.refs[nm][off:off+ln]
        # The simulated read can't overlap a non-A-C-G-T character in the
        # reference
        while self.__re.match(seq):
            off = random.randint(0, self.lens[refi] - ln)
            seq = self.refs[nm][off:off+ln]
        # Pick some reads to reverse-complement
        if random.random() > 0.5:
            fw = False
            seq = revcomp(seq)
        if verbose:
            print >>sys.stderr, "...done"
        return (nm, off, fw, seq)

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
    def __init__(self, inp, rddist1, rddist2):
        self.inp = inp
        self.rddist1, self.rddist2 = rddist1, rddist2
    
    def __iter__(self):
        for (rd1, rd2) in self.inp:
            assert rd1 is not None
            self.rddist1.addRead(rd1)
            if rd2 is not None:
                self.rddist2.addRead(rd2)
            #
            yield (rd1, rd2)

class Output(threading.Thread):
    """ Encapsulates the output reader.  Reads SAM output from the aligner,
        updates empirical distributions for e.g. fragment length, and looks for
        records that correspond to simulated reads.  If a record corresponds to
        a simulated read, its correctness will be checked and a tuple
        written """
    
    def __init__(self, samIfh, samOfh, trainFh, scoreDist, fragDist):
        threading.Thread.__init__(self)
        self.samIfh = samIfh   # SAM records come from here
        self.samOfh = samOfh   # normal (non-simulated) SAM records go here
        self.trainFh = trainFh # write training data here
        self.fragDist = fragDist
        self.scoreDist = scoreDist
    
    def run(self):
        lastAl = None
        for ln in self.samIfh:
            if ln.startswith('@'):
                self.samOfh.write(ln) # header line
                continue
            al = Alignment(ln)
            nm, flags, refid, pos, _, _, _, _, _, seq, qual, _ = string.split(ln, '\t', 11)
            if al.name.startswith('!!ts!!'):
                # this is a simulated read
                _, chr, fw, off = string.split(al.name, '_')
                off = int(off)
                if al.isAligned():
                    correct = False
                    # Check reference id, orientation
                    if chr == al.refid and fw == al.orientation():
                        # Check offset
                        correct = abs(off - al.pos) < args.wiggle
                    # Output training tuple:
                    # 1. Read length
                    # 2. Alignment type (UU, UP, CP, DP)
                    # 3. Score of best
                    # 4. Score of second-best
                    # 5. Fragment length (if paired-end)
                    # 6. Mate's score of best
                    # 7. Mate's score of second-best
                    # 8. Correct?
                    trainFh.write("\t".join([str(len(al)), al.alType,
                                  str(al.bestScore), str(al.secondBestScore),
                                  str(int(correct))]))
            else:
                # Take alignment info into account
                if lastAl is not None and al.name == lastAl.name:
                    mate1, mate2 = al, lastAl
                    if (lastAl.flags & 64) != 0:
                        mate1, mate2 = mate2, mate1
                    alPair = AlignmentPair(mate1, mate2)
                    self.fragDist.addFragment(alPair) # alignment pair
                if al.isAligned():
                    self.scoreDist.addAlignment(al)
                # Send SAM to SAM output filehandle
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
                self.fh.close()
                break
            self.fh.write(item)
            self.fh.flush()
            self.q.task_done()
            i += 1

class Aligner(object):
    """ Encapsulates the aligner """
    def __init__(self):
        pass

def createInput(rddist1, rddist2):
    """ Return an Input object that reads all user-provided input reads """
    if args.fastq:
        return InputWrapper(Input(format="fastq", unpFns=args.U, m1Fns=args.m1, m2Fns=args.m2), rddist1, rddist2)
    elif args.fasta:
        return InputWrapper(Input(format="fasta", unpFns=args.U, m1Fns=args.m1, m2Fns=args.m2), rddist1, rddist2)

def go():
    """ Main driver for tandem simulator """
    bt2_cmd = "bowtie2 "
    if args.bt2_exe is not None:
        bt2_cmd = args.bt2_exe + " "
    if args.bt2_args is not None:
        bt2_cmd += args.bt2_args
    fifoArgs = ""
    fifoFns = [None, None, None]
    fifoFhs = [None, None, None]
    fifoQs  = [None, None, None]
    fifoWs  = [None, None, None]
    tmpdir = tempfile.mkdtemp()
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
    bt2_cmd += " --reorder --sam-no-qname-trunc -q " + fifoArgs
    if args.verbose:
        print >> sys.stderr, "Bowtie 2 command: '%s'" % bt2_cmd
    pipe = subprocess.Popen(bt2_cmd, shell=True, stdout=subprocess.PIPE)
    if args.U  is not None: fifoFhs[0] = open(fifoFns[0], 'w')
    if args.m1 is not None: fifoFhs[1] = open(fifoFns[1], 'w')
    if args.m2 is not None: fifoFhs[2] = open(fifoFns[2], 'w')
    scDist = ScoreDist()
    fragDist = FragmentLengthDist()
    samOfh = sys.stdout
    if args.S:
        samOfh = open(args.S, 'w')
    othread = Output(pipe.stdout, samOfh, sys.stderr, scDist, fragDist)
    othread.start()
    rddist1, rddist2 = ReadLengthDist(), ReadLengthDist()
    sctok = map(int, string.split(args.scoring, ','))
    scoring = Scoring(sctok[0], (sctok[1], sctok[2]), sctok[3], (sctok[4], sctok[5]), (sctok[6], sctok[7]))
    for i in xrange(0, 3):
        if fifoFhs[i] is not None:
            fifoQs[i] = Queue()
            fifoWs[i] = AsyncWriter(fifoFhs[i], fifoQs[i], "Thread %d" % i)
            fifoWs[i].start()
    for (rd1, rd2) in iter(createInput(rddist1, rddist2)):
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
    for i in xrange(0, 3):
        if fifoFhs[i] is not None:
            fifoQs[i].put(None)
    for i in xrange(0, 3):
        if fifoFhs[i] is not None:
            fifoQs[i].join()
            fifoWs[i].join()
            fifoFhs[i].close()
            os.unlink(fifoFns[i])
    othread.join()
    if args.S: samOfh.close()
    print "rddist1:"
    print rddist1
    print "rddist2:"
    print rddist2
    print "scDist:"
    print scDist
    print "fragDist:"
    print fragDist

if __name__ == "__main__":

    if args.test:
        import unittest
        class Test(unittest.TestCase):
            def test_1(self):
                pass
        unittest.main(argv=[sys.argv[0]])
    elif args.profile:
        import cProfile
        cProfile.run('go()')
    else:
        go()
