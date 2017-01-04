#!/bin/sh

awk -v OFS=',' '{print $0, FILENAME}' *.csv > combined.csv
