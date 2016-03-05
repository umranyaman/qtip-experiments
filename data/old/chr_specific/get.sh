#!/bin/sh

mkdir -p hs_chr21 hs_chr15 hs_X_exome mm_17
python ../../bin/downloader.py --manifest=MANIFEST
