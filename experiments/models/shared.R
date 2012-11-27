library(bitops)
library(verification)

# Open a tabulated sam file
openTabulatedSam <- function(fn) {
	return(read.table(fn, comment.char="", quote="", header=T))
}

# Given data frame of alignments, return data frame with just aligned reads
selectAligned <- function(x) {
	return(bitAnd(x$flag, 4) == 0)
}

# Given data frame of alignments, return just those that aligned once
# (i.e. AS.i is defined but XS.i is not)
selectAlignedOnce <- function(x) {
	return(bitAnd(x$flag, 4) == 0 & is.na(x$XS.i))
}

# Given data frame of alignments, return just those that aligned twice
# or more (i.e. AS.i and XS.i are defined)
selectAlignedMoreThanOnce <- function(x) {
	return(bitAnd(x$flag, 4) == 0 & !is.na(x$XS.i))
}

# Standardize the given vector so that all the elements are on
# [0.0, 1.0]
bt0and1 <- function(x) {
	x <- x - min(x[!is.na(x)])
	x <- x / max(x[!is.na(x)])
	return(x)
}

# Assess quality of ranking.  In this case we add up all absolute
# differences between element at rank i, and N/i where N is the total
# number of elements.
rankingError1 <- function(x) {
	ordr <- order(x$mapq)
	dif <- x$ZC.i[ordr] - seq(0, 1, length.out=nrow(x))
	return(sum(abs(dif)))
}

# Another way to assess the quality of ranking.  Here for each
# incorrect alignment we add its rank to the penalty, where the
# worst-ranked alignment has rank 1.
rankingError2 <- function(x) {
	ordr <- order(x$mapq)
	return(sum(which(x$ZC.i[ordr] == 0)))
}

# Plot a ROC.  Not really appropriate for comparing two different sets
# of alignments, since the "hit rate" and "false alarm rate" ignore the
# fact that one set of alignments might contain more/different
# alignments than the other.  A similar function that simply plots the
# number of correct alignments on the vertical axis and number of
# incorrect alignments on the horizontal axis would be better. 
plotRoc <- function(x) {
	roc.plot(x$ZC.i, x$mapq)
}
