source("shared.R")

# Given data frame containing alignments that aligned just once, return
# a fit for a linear model with response = mapq ranking and terms =
# alignment score.
modelAlignedOnce <- function(x) {
	model <- lm(ZC.i ~ AS.i, x)
	mapq <- model$coefficients[1] + model$coefficients[2] * x$AS.i
	return(list(model=model, mapq=mapq))
}

# Given data frame containing alignments that aligned just once, return
# a fit for a linear model with response = mapq ranking and terms =
# alignment score.
modelAlignedMoreThanOnce <- function(x) {
	x$AS.iMXS.i <- x$AS.i - x$XS.i
	model <- lm(ZC.i ~ AS.i + AS.iMXS.i, x)
	mapq <- model$coefficients[1] + model$coefficients[2] * x$AS.i + model$coefficients[3] * x$AS.iMXS.i
	return(list(model=model, mapq=mapq))
}

# Given a data frame containing all aligned reads (including those that
# align just once or many times) return a fit for a linear model with
# response = mapq ranking and terms = alignment scores.
#
# Second argument is a scaling factor, which determines how XS.i values
# should be assigned to reads that align only once.  Right now we set
# them equal to max(AS.i) - 2 * (max(AS.i) - min(AS.i)).  
modelAltogether <- function(x, sca=2.0) {
	repl <- max(x$AS.i) - sca * (max(x$AS.i) - min(x$AS.i))
	xs <- ifelse(is.na(x$XS.i), repl, x$XS.i)
	x$AS.iMXS.i <- x$AS.i - xs
	model <- lm(ZC.i ~ AS.i + AS.iMXS.i, x)
	mapq <- model$coefficients[1] + model$coefficients[2] * x$AS.i + model$coefficients[3] * x$AS.iMXS.i
	return(list(model=model, mapq=mapq))
}

# Plot the ZC.i values sorted by mapq, along with mapqs
plotLinesAndDots <- function(x, plotReps=F) {
	ordr <- order(x$mapq)
	plot(jitter(x$ZC.i[ordr]), col=rgb(0, 0, 0, 0.1))
	mapq <- bt0and1(x$mapq[ordr])
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

# Plot a histogram of the mapqs for the incorrectly aligned reads. 
incorrectMapqHist <- function(x) { hist(x$mapq[x$ZC.i == 0]) }

# Return a table of the mapqs for the incorrectly aligned reads 
incorrectMapqTable <- function(x) { return(table(x$mapq[x$ZC.i == 0])) }

# Given filenames for SAM files emitted by two tools, analyze 
fitMapqModels <- function(x) {
	tab <- openTabulatedSam(x)
	tab.1 <- tab[selectAlignedOnce(tab),]
	tab.2 <- tab[selectAlignedMoreThanOnce(tab),]
	tab.all <- rbind(tab.1, tab.2)
	fit.1 <- modelAlignedOnce(tab.1); tab.1$mapq <- bt0and1(fit.1$mapq)
	fit.2 <- modelAlignedMoreThanOnce(tab.2); tab.2$mapq <- bt0and1(fit.2$mapq)
	fit.all <- modelAltogether(tab.all, sca=2.0); tab.all$mapq <- bt0and1(fit.all$mapq)
	rankingError(tab.all)
	rankingError(rbind(tab.1, tab.2))
}

x <- "/tmp/r0_ill_100_100k.bt2_s.sat"
