source("shared.R")

#
# Problems with this model:
#
# 1. It's clear from plotLinesAndDots(tab.all, replace.na=T) that the
#    unique alignments mostly segregate to the right (higher mapqs)
#    which causes some of the incorrectly-aligned unique alignments to
#    receive high ranks.  At the same time, fiddling with the 'sca'
#    parameter to plotLinesAndDots (which determines how we replace the
#    "NA" second-best alignment scores) doesn't do much to reduce
#    ranking error.
# 2. If we fit separate models for the unique alignments and for the
#    non-unique alignments, we might get more sensible rankings within
#    each subset, but then we don't know how to combine them into one
#    overall ranking, which is needed to compare models.
# 3. It's clear from plotLinesAndDots(tab.2) that the alignments where
#    AS.i - XS.i = 0 are in a special class, and using
#    modelAlignedMoreThanOnce with exp=1 seems to rank them too high
#    (they appear together in a chunk in the middle/left).
#    Experimenting with setting exp by hand did seem to improve this.
#    exp=0.8 mostly shifted the AS.i - XS.i = 0 alignments to the left
#    edge and improved ranking error.
#

# Given a data frame containing all aligned reads (including those that
# align just once or many times) return a fit for a linear model with
# response = mapq ranking and terms = alignment scores.
#
# Second argument is a scaling factor, which determines how XS.i values
# should be assigned to reads that align only once.  Right now we set
# them equal to max(AS.i) - 2 * (max(AS.i) - min(AS.i)).  
modelAllLogistic <- function(x, replace.na=F, sca=2.0) {
	if(replace.na) {
		repl <- max(x$AS.i) - sca * (max(x$AS.i) - min(x$AS.i))
		xs <- ifelse(is.na(x$XS.i), repl, x$XS.i)
	} else {
		xs <- x$XS.i
	}
	x$AS.iMXS.i <- x$AS.i - xs # could have some NAs
	model <- glm(ZC.i ~ AS.i + AS.iMXS.i, data=x, family=binomial("logit"))
	mapq <- model$coefficients[1] + model$coefficients[2] * x$AS.i + model$coefficients[3] * x$AS.iMXS.i
	return(list(model=model, mapq=mapq))
}

# Given filenames for SAM files emitted by two tools, analyze 
fitMapqModelsLogistic <- function(x) {
	tab <- openTabulatedSam(x)
	tab.all <- tab[selectAligned(tab),]
	fit.all <- modelAllLogistic(tab.all, replace.na=T, sca=1.2)
	tab.all$model_mapq <- bt0and1(fit.all$mapq)
	return(rankingError2(tab.all))
}

if(False) {
	dr <- "/Users/langmead/Documents/workspace/mapq"
	x <- paste(dr, "examples/mason/r0_ill_100_100k.bt2_s.sat.bz2", sep="/")
	fitMapqModelsLogistic(x)
}
