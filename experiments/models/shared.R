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

# Another way to assess the quality of ranking.  Here for each
# incorrect alignment we add its rank to the penalty, where the
# worst-ranked alignment has rank 1.
rankingError <- function(x, mapq) {
	ordr <- order(mapq)
	return(sum(which(x$ZC.i[ordr] == 0)))
}

# Another way to assess the quality of ranking.  Here for each
# correct alignment we add its rank to the score.  Higher is better.
rankingQuality1 <- function(x) {
	ordr <- order(x$mapq)
	ordr_model <- order(x$model_mapq)
	return(c(sum(which(x$ZC.i[ordr] == 1)/length(x$ZC.i)), sum(which(x$ZC.i[ordr_model] == 1)/length(x$ZC.i))))
}

# Plot the ZC.i values sorted by mapq, along with mapqs
plotLinesAndDots <- function(x, mapq, plotReps=F) {
	ordr <- order(mapq)
	plot(jitter(x$ZC.i[ordr]), col=rgb(0, 0, 0, 0.1))
	mapq <- bt0and1(mapq[ordr])
	points(mapq, col=rgb(1.0, 0.0, 0.0, 0.1))
	xsi <- bt0and1(x$XS.i[ordr])
	asi <- bt0and1(x$AS.i[ordr])
	points(asi - xsi, col=rgb(0.0, 0.0, 1.0, 0.1))
	points(asi, col=rgb(1.0, 0.5, 0.25, 0.1))
	if(plotReps) {
		re <- bt0and1(x$AllRepeats[ordr])
		points(re, col=rgb(0.0, 0.5, 0.0, 0.1))
	}
}

# Among the incorrect alignments, return the top n according to
# highest mapping quality.
topNincorrect <- function(x, mapq, n) {
	x.incor <- x[x$ZC.i == 0,]
	mapq.incor <- mapq[x$ZC.i == 0]
	ordr <- order(mapq.incor, decreasing=T)
	return(head(x.incor[ordr,], n=n))
}

# Plot a histogram of the mapqs for the incorrectly aligned reads. 
incorrectMapqHist <- function(x) { hist(x$model_mapq[x$ZC.i == 0]) }

# Return a table of the mapqs for the incorrectly aligned reads 
incorrectMapqTable <- function(x) { return(table(x$model_mapq[x$ZC.i == 0])) }
