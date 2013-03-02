#!/usr/bin/perl

use strict;
use warnings;

# Given FASTQ output from mason -i, generate a FASTQ file with Art-style
# read names.

my $fq_fn = shift @ARGV;

open(my $fq_fh, $fq_fn)  || die "Could not open '$fq_fn' for reading";

while(my $name  = readline $fq_fh) {
	my $seq   = readline $fq_fh;
	my $name2 = readline $fq_fh;
	my $qual  = readline $fq_fh;
	defined($qual) || die "Files ended at different reads";
	$name =~ /orig_begin=([0-9]*)/;
	my $left = $1;
	$name =~ /orig_end=([0-9]*)/;
	my $right = $1;
	$name =~ /contig=([^\s]*)/;
	my $chr = $1;
	$left < $right || die;
	my $fw = ($name =~ /strand=forward/) ? "+" : "-";
	my $newname = "\@$chr:$left:$fw\n";
	print "$newname$seq$name2$qual";
}

close($fq_fh);
