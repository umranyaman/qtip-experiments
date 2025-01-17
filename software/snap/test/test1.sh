#!/bin/sh


#  -om  Output multiple alignments.  Takes as a parameter the maximum extra edit distance relative to the best alignment
#       to allow for secondary alignments
# -omax Limit the number of alignments per read generated by -om.  This means that if -om would generate more
#       than -omax secondary alignments, SNAP will write out only the best -omax of them, where 'best' means
#       'with the lowest edit distance'.  Ties are broken arbitrarily.

# For "paired" mode:
#   -s   min and max spacing to allow between paired ends (default: 50 1000).
#   -fs  force spacing to lie between min and max.

BASE_URL=http://www.cs.jhu.edu/~langmea/resources

[ ! -f lambda_virus.fa ] && wget $BASE_URL/lambda_virus.fa
[ ! -f reads_1.fq ] && wget $BASE_URL/reads_1.fq
[ ! -f reads_2.fq ] && wget $BASE_URL/reads_2.fq

# Index
../snap/snap-aligner index lambda_virus.fa lambda_virus.snap

# Simple unpaired
../snap/snap-aligner single lambda_virus.snap reads_1.fq -o -sam unp1.sam

# Simple unpaired with omax
../snap/snap-aligner single lambda_virus.snap reads_1.fq -o -sam unp1_omax1.sam -om -omax 1

# Unpaired using stdin
cat reads_1.fq | ../snap/snap-aligner single lambda_virus.snap -fastq - -o -sam unp2.sam

# Paired-end, separate files
../snap/snap-aligner paired lambda_virus.snap reads_1.fq reads_2.fq -o -sam pai1.sam

python ../../../bin/fastq_interleave.py reads_1.fq reads_2.fq > reads_paired.fq

# Paired-end, interleaved file
../snap/snap-aligner paired lambda_virus.snap -pairedInterleavedFastq reads_paired.fq -o -sam pai2.sam

# Paired-end, interleaved, stdin file
cat reads_paired.fq | ../snap/snap-aligner paired lambda_virus.snap -pairedInterleavedFastq - -o -sam pai3.sam

# Paired-end twice
../snap/snap-aligner paired lambda_virus.snap reads_1.fq reads_2.fq reads_1.fq reads_2.fq -o -sam pai_d1.sam

# Paired-end fragment length experiments

cat <<EOF >.r1.fq
@r1
GGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

cat <<EOF >.r2.fq
@r1
GCCTCGCTTTCAGCACCTGTCGTTTCCTTTCTTTTCAGAGGGTATTTTAAATAAAAACATTAAGTTATGA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

# Just some random bases
cat <<EOF >.r2.bad.fq
@r1
ACGTATTATATGGCCGTAGAGTGCGCATAGATGCTCAGTCAAACCCGCGGATATATAAAACCGGCGGTTT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

# These are the circumstances that yield a concordant alignment
../snap/snap-aligner paired lambda_virus.snap .r1.fq .r2.fq -s 69 70 -o -sam paired_conc_1.sam

# This gives a discordant alignment
../snap/snap-aligner paired lambda_virus.snap .r1.fq .r2.fq -s 1 2 -fs -o -sam paired_disc_1.sam

# This gives a pair where just one end aligns
../snap/snap-aligner paired lambda_virus.snap .r1.fq .r2.bad.fq -s 69 70 -o -sam paired_unp_1.sam

# Now let's try to figure out exactly what distances the -s parameter is setting.
# Theory: it's specifying a right-open interval where the interval describes how
# far the leftmost aligned base of the second mate is allowed to be from the
# leftmost aligned base of the first mate.

cat <<EOF >.r2.nudged_right.fq
@r1
GCCTCGCTTTCAGCACCTGTCGTTTCCTTTCTTTTCAGAGGGTATTTTAAATAAAAACATTAAGTTATG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

../snap/snap-aligner paired lambda_virus.snap .r1.fq .r2.nudged_right.fq -s 70 71 -fs -o -sam paired_conc_2.sam

cat <<EOF >.r2.nudged_left.fq
@r1
GCCTCGCTTTCAGCACCTGTCGTTTCCTTTCTTTTCAGAGGGTATTTTAAATAAAAACATTAAGTTATGAC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

# yep, this gives concordant alignment
../snap/snap-aligner paired lambda_virus.snap .r1.fq .r2.nudged_left.fq -s 68 69 -fs -o -sam paired_conc_3.sam
