"""
tandem_sam_scores.py

For each alignment, compare the "target" simulated alignment score to the
actual score obtained by the aligner.  When the read is simulated, we borrow
the target score and the pattern of mismatches and gaps from an input
alignment.  But because the new read's sequence and point of origin are
different, and because the aligner's heuristics might act differently on the
tandem read than on the input read, the aligned score might be different.
Could be either higher or lower.  Here we compare and make a table showing
how scores change before and after.
"""

from __future__ import print_function
import sys
from collections import defaultdict

scores = defaultdict(int)

for ln in sys.stdin:
    if ln[0] == '@':
        continue
    toks = ln.rstrip().split('\t')
    assert toks[0][:6] == '!!ts!!'
    ref_id, fw, ref_off, expected_score, typ = toks[0][16:].split('!!ts-sep!!')
    actual_score = None
    for ef in toks[11:]:
        if ef.startswith('AS:i:'):
            actual_score = int(ef[5:])
            break
    scores[(expected_score, actual_score)] += 1

for k, v in sorted(scores.items()):
    expected_score, actual_score = k
    actual_score = 'NA' if actual_score is None else str(actual_score)
    print("%s,%s,%d" % (expected_score, actual_score, v))
