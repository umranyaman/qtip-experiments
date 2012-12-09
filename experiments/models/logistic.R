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
modelAllLogistic <- function(x, do.rpart=F, incl.Xs=F, incl.repeats=F, incl.xd=F, replace.na=F, sca=2.0) {
	if(replace.na) {
		repl <- max(x$AS.i) - sca * (max(x$AS.i) - min(x$AS.i))
		xs <- ifelse(is.na(x$XS.i), repl, x$XS.i)
		if(incl.Xs) {
			xs <- pmax(xs, x$Xs.i, na.rm=T)
		}
	} else {
		xs <- x$XS.i
	}
	x$AS.iMXS.i <- x$AS.i - xs # could have some NAs
	form <- "x$ZC.i ~ x$AS.i + x$AS.iMXS.i"
	if(incl.repeats) {
		form <- paste(form, "+ x$AllRepeats")
	}
	if(incl.xd) {
		form <- paste(form, "+ x$XD.i")
	}
	form <- as.formula(form)
	if(do.rpart) {
		require('rpart')
		model <- rpart(form)
	} else {
		model <- glm(form, family=binomial("logit"))
	}
	mapq <- model$coefficients[1] + model$coefficients[2] * x$AS.i + model$coefficients[3] * x$AS.iMXS.i
	return(list(model=model, mapq=mapq))
}

bestSca <- function(tab, do.rpart=F, incl.Xs=F, incl.repeats=F, incl.xd=F, sca=1.2, maxit=25, method="SANN") {
	fn <- function(x) {
		fit <- modelAllLogistic(
			tab,
			do.rpart=do.rpart,
			incl.Xs=incl.Xs,
			incl.repeats=incl.repeats,
			incl.xd=incl.xd,
			replace.na=T,
			sca=x)
		err <- rankingError(tab, fit$mapq)
		print(paste("Called with x=", x, " error=", err, sep=""))
		return(err)
	}
	opt <- optim(sca, fn, method=method, lower=0.5, upper=5.0, control=list(maxit=maxit, reltol=0.01, trace=T))
	return(opt)
}

# Given filenames for SAM files emitted by two tools, analyze 
fitMapqModelsLogistic <- function(x, do.rpart=F, incl.Xs=T, incl.repeats=F, incl.xd=F, sca=1.2) {
	tab <- openTabulatedSam(x)
	tab.all <- tab[selectAligned(tab),]
	opt <- bestSca(tab.all, method="Nelder-Mead", maxit=1, sca=sca)
	sca <- opt$par
	fit.all <- modelAllLogistic(tab.all, do.rpart=do.rpart, incl.Xs=incl.Xs, incl.repeats=incl.repeats, incl.xd=incl.xd, replace.na=T, sca=sca)
	tab.all$model_mapq <- bt0and1(fit.all$mapq)
	return(list(tb=tab.all, fit=fit.all, err1=rankingError(tab.all, tab.all$model_mapq), err2=rankingError(tab.all, tab.all$mapq)))
}

if(F) {
	dr <- "/Users/langmead/Documents/workspace/mapq"
	x <- paste(dr, "examples/mason/r0_ill_100_100k.bt2_s.sat.bz2", sep="/")
	fit <- fitMapqModelsLogistic(x)
	roc_table_compare(fit$tb, fit$tb$model_mapq, fit$tb$mapq)
}
