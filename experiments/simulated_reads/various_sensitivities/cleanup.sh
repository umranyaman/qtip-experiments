#!/bin/sh

# For finding big files:
# find . -type f -printf "%s\t%p\n" | sort -n | tail -n 60

# remove temporaries
rm -rf *.temp

# remove all but one experiment
#find . -name 'input_intermediates_rec_*.npy' | grep -v '05m_30s' | xargs rm -f

# remove all but one trial
find . -name 'final.sam' | grep -v 'trial0' | xargs rm -f
find . -name 'tandem_*.sam' | grep -v 'trial0' | xargs rm -f
find . -name 'predictions.*.npy' | grep -v 'trial0' | xargs rm -f
find . -name 'predictions_assess.*.npy' | grep -v 'trial0' | xargs rm -f
find . -name 'tandem_intermediates_reads_*.fastq' | grep -v 'trial0' | xargs rm -f
find . -name 'tandem_intermediates_rec_*.npy' | grep -v 'trial0' | xargs rm -f
