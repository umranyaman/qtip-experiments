#!/usr/bin/env python

"""
table_f.py

Scan all the ROCs output by vcfroc and record, for two different levels of
stratification, the best f-beta-score in that stratum and the cutoffs/other
stats associated with that best. Do this for various values of beta in the
f-beta score.
"""

from __future__ import print_function


def f_beta(precision, recall, beta):
    if precision + recall == 0:
        return 0
    return (1 + beta * beta) * precision * recall / (beta * beta * precision + recall)


def go():
    nm = 'ERR194147'
    betas = list(map(lambda x: 2 ** (x/10.0), range(-10, 11)))
    betas = list(sorted(betas + [0.75, 1.5]))
    beta_labs = list(map(lambda x: ('%0.3f' % x).replace('.', '_'), betas))

    # Print headers
    print(','.join(['name', 'filt', 'chrom', 'regime', 'mapq_cutoff']), end='')
    for lab in beta_labs:
        print(',' + ','.join(['f_' + lab,
                              'qual_cutoff_' + lab,
                              'mapq_cutoff_' + lab,
                              'regime_' + lab,
                              'precision_' + lab,
                              'recall_' + lab,
                              'tp_' + lab,
                              'fp_' + lab,
                              'fn_' + lab]), end='')
    print()

    for ch in [1, 6, 19, 22]:
        for filt in ['cr_filt', 'rmsk_filt']:

            best_mapq = [0] * len(betas)
            argbest_mapq = [None] * len(betas)

            for regime in ['input', 'final']:

                best_mapq_r = [0] * len(betas)
                argbest_mapq_r = [None] * len(betas)

                for mapq_cutoff in ['00', '01', '02', '03', '05', '08', '10', '15', '20', '30', 'd', 's', 'u']:

                    best_mapq_rc = [0] * len(betas)
                    argbest_mapq_rc = [None] * len(betas)

                    fn = '%s.sam/%s_%s_%d_%s.%s.roc' % (nm, nm, regime, ch, mapq_cutoff, filt)
                    with open(fn) as fh:
                        for ln in fh:
                            toks = ln.rstrip().split('\t')
                            if toks[0] == 'threshold':
                                continue
                            threshold = float(toks[0])
                            num_snps, false_positive_snps, false_negative_snps = map(int, toks[1:4])
                            tp, fp, fn = num_snps - false_positive_snps, false_positive_snps, false_negative_snps
                            precision, recall = 0, 0
                            if tp + fp > 0:
                                precision = float(tp) / (tp+fp)
                            if tp + fn > 0:
                                recall = float(tp) / (tp+fn)
                            f_scores = list(map(lambda b: f_beta(precision, recall, b), betas))
                            for i, f in enumerate(f_scores):
                                # Better than best so far for this chromosome, MAPQ regime and MAPQ cutoff?
                                if f >= best_mapq_rc[i]:
                                    best_mapq_rc[i] = f
                                    argbest_mapq_rc[i] = (threshold, mapq_cutoff, regime, precision, recall, tp, fp, fn)
                                # Better than best so far for this chromosome and MAPQ regime?
                                if f >= best_mapq_r[i]:
                                    best_mapq_r[i] = f
                                    argbest_mapq_r[i] = (threshold, mapq_cutoff, regime, precision, recall, tp, fp, fn)
                                # Better than best so far for this chromosome?
                                if f >= best_mapq[i]:
                                    best_mapq[i] = f
                                    argbest_mapq[i] = (threshold, mapq_cutoff, regime, precision, recall, tp, fp, fn)

                    print(','.join(map(str, [nm, filt, ch, regime, mapq_cutoff])), end='')
                    for f, tup in zip(best_mapq_rc, argbest_mapq_rc):
                        qual_cutoff, mapq_cutoff, regime, precision, recall, tp, fp, fn = tup
                        print(',' + ','.join(map(str, [f, qual_cutoff, mapq_cutoff, regime, precision, recall, tp, fp, fn])), end='')
                    print()

                print(','.join(map(str, [nm, filt, ch, regime, '*'])), end='')
                for f, tup in zip(best_mapq_r, argbest_mapq_r):
                    qual_cutoff, mapq_cutoff, regime, precision, recall, tp, fp, fn = tup
                    print(',' + ','.join(map(str, [f, qual_cutoff, mapq_cutoff, regime, precision, recall, tp, fp, fn])), end='')
                print()

            print(','.join(map(str, [nm, filt, ch, '*', '*'])), end='')
            for f, tup in zip(best_mapq, argbest_mapq):
                qual_cutoff, mapq_cutoff, regime, precision, recall, tp, fp, fn = tup
                print(',' + ','.join(map(str, [f, qual_cutoff, mapq_cutoff, regime, precision, recall, tp, fp, fn])), end='')
            print()

if __name__ == "__main__":
    go()
