#!/usr/bin/env python

"""
global_align.py

Functions for global alignment, with and without affine gap penalties.
"""

__author__ = "Ben Langmead"
__email__ = "langmea@cs.jhu.edu"

import numpy

def globalAlign(s1, s2, subs, s1g, s2g, backtrace=False):
    
    """ Perform Needleman-Wunsch-style global alignment of s1 (rows) and s2
        (columns). """
    
    nrow, ncol = len(s1), len(s2)
    tab = numpy.zeros((nrow+1, ncol+1), int)
    bt_dg, bt_up, bt_lf = None, None, None
    if backtrace:
        bt_dg = numpy.zeros((nrow, ncol), dtype=bool)
        bt_up = numpy.zeros((nrow, ncol), dtype=bool)
        bt_lf = numpy.zeros((nrow, ncol), dtype=bool)
    for i in xrange(1, nrow+1): tab[i][0] = i * s2g
    for j in xrange(1, ncol+1): tab[0][j] = j * s1g
    for i in xrange(1, nrow+1):
        for j in xrange(1, ncol+1):
            up = tab[i-1, j  ] + s2g
            lf = tab[i  , j-1] + s1g
            dg = tab[i-1, j-1] + subs(s1[i-1], s2[j-1])
            tab[i, j] = min(up, lf, dg)
            if backtrace:
                bt_up[i-1, j-1] = up == tab[i, j]
                bt_lf[i-1, j-1] = lf == tab[i, j]
                bt_dg[i-1, j-1] = dg == tab[i, j]
    return (tab[len(s1)][len(s2)], tab, (bt_dg, bt_up, bt_lf))

def globalAlignAffine(s1, s2, subs, s1g_op, s1g_ex, s2g_op, s2g_ex, backtrace=False):
    
    """ Perform Needleman-Wunsch-style global alignment of s1 (rows) and s2
        (columns).  Use affine gap penalty. """
    
    nrow, ncol = len(s1), len(s2)
    E, F, H = numpy.zeros((nrow+1, ncol+1), int), \
              numpy.zeros((nrow+1, ncol+1), int), \
              numpy.zeros((nrow+1, ncol+1), int)
    bt_h_dg, bt_h_up_op, bt_h_up_ex, bt_h_lf_op, bt_h_lf_ex = \
        None, None, None, None, None
    bt_e_lf_op, bt_e_lf_ex = None, None
    bt_f_up_op, bt_f_up_ex = None, None
    if backtrace:
        bt_h_dg    = numpy.zeros((nrow, ncol), dtype=bool)
        bt_h_up_op = numpy.zeros((nrow, ncol), dtype=bool)
        bt_h_up_ex = numpy.zeros((nrow, ncol), dtype=bool)
        bt_h_lf_op = numpy.zeros((nrow, ncol), dtype=bool)
        bt_h_lf_ex = numpy.zeros((nrow, ncol), dtype=bool)
        bt_e_lf_op = numpy.zeros((nrow, ncol), dtype=bool)
        bt_e_lf_ex = numpy.zeros((nrow, ncol), dtype=bool)
        bt_f_up_op = numpy.zeros((nrow, ncol), dtype=bool)
        bt_f_up_ex = numpy.zeros((nrow, ncol), dtype=bool)
    for i in xrange(1, nrow+1):
        H[i][0] = F[i][0] = s2g_op + (i * s2g_ex)
    for j in xrange(1, ncol+1):
        H[0][j] = E[0][j] = s1g_op + (j * s1g_ex)
    for i in xrange(1, nrow+1):
        for j in xrange(1, ncol+1):
            h_up_op = H[i-1, j  ] + s2g_op + s2g_ex
            h_up_ex = F[i-1, j  ] + s2g_ex
            h_lf_op = H[i  , j-1] + s1g_op + s1g_ex
            h_lf_ex = E[i  , j-1] + s1g_ex
            #h_dg    = H[i-1, j-1] + subs(s1[i-1], s2[j-1])
            h_dg    = H[i-1, j-1] + (0 if s1[i-1] == s2[j-1] else 6)
            f_up_op = H[i-1, j  ] + s2g_op + s2g_ex
            f_up_ex = F[i-1, j  ] + s2g_ex
            e_lf_op = H[i  , j-1] + s1g_op + s1g_ex
            e_lf_ex = E[i  , j-1] + s1g_ex
            H[i, j] = min(h_up_op, h_up_ex, h_lf_op, h_lf_ex, h_dg)
            E[i, j] = min(e_lf_op, e_lf_ex)
            F[i, j] = min(f_up_op, f_up_ex)
            if backtrace:
                bt_h_dg   [i-1, j-1] = h_dg    == H[i, j]
                bt_h_up_op[i-1, j-1] = h_up_op == H[i, j]
                bt_h_up_ex[i-1, j-1] = h_up_ex == H[i, j]
                bt_h_lf_op[i-1, j-1] = h_lf_op == H[i, j]
                bt_h_lf_ex[i-1, j-1] = h_lf_ex == H[i, j]
                bt_e_lf_op[i-1, j-1] = e_lf_op == E[i, j]
                bt_e_lf_ex[i-1, j-1] = e_lf_ex == E[i, j]
                bt_f_up_op[i-1, j-1] = f_up_op == F[i, j]
                bt_f_up_ex[i-1, j-1] = f_up_ex == F[i, j]
    return (H[len(s1)][len(s2)], (H, E, F), (bt_h_dg, bt_h_up_op, bt_h_up_ex, bt_h_lf_op, bt_h_lf_ex, bt_e_lf_op, bt_e_lf_ex, bt_f_up_op, bt_f_up_ex))
