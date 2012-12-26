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
# align just once or many times) return a fit for a logistic regression
# model with response = mapq ranking and terms = alignment scores.
#
# "sca" is a scaling factor, which determines how XS.i values should be
# assigned to reads that align only once.  Right now we set them equal
# to max(AS.i) - 2 * (max(AS.i) - min(AS.i)).  
modelAllLogistic <- function(x, incl.Xs=F, incl.repeats=F, incl.xd=F, incl.len=F, replace.na=T, rescale=T, sca=2.0) {
	cov <- getCovars(x, incl.Xs=incl.Xs, incl.xd=incl.xd, replace.na=replace.na, rescale=rescale, sca=sca)
	form <- "x$ZC.i ~ cov$as + cov$xs"
	if(incl.repeats) { form <- paste(form, "+ x$AllRepeats") }
	if(incl.xd)      { form <- paste(form, "+ x$XD.i") }
	if(incl.len)     { form <- paste(form, "+ nchar(x$seq)") }
	model <- glm(as.formula(form), family=binomial("logit"))
	coef <- coefficients(model)
	mapq <- 1.0 / (1.0 + exp(-(coef[1] + coef[2] * cov$as + coef[3] * cov$xs)))
	mapq2 <- coef[1] + coef[2] * cov$as + coef[3] * cov$xs
	return(list(model=model, mapq=mapq, mapq2=mapq2, coef=coef))
}

# Use "optim" to search for the best scaling factor to use for the dataset
bestSca <- function(tab, incl.Xs=F, incl.repeats=F, incl.xd=F, incl.len=F, rescale=T, sca=1.2, maxit=60, method="SANN") {
	fn <- function(x) {
		fit <- modelAllLogistic(
			tab,
			incl.Xs=incl.Xs,
			incl.repeats=incl.repeats,
			incl.xd=incl.xd,
			incl.len=incl.len,
			replace.na=T,
			rescale=rescale,
			sca=x)
		err <- rankingError(tab, fit$mapq)
		print(paste("Called with x=", x, " error=", err, sep=""))
		return(err)
	}
	opt <- optim(sca, fn, method=method, control=list(maxit=maxit))
	return(opt)
}

# Given filenames for SAM files emitted by two tools, analyze 
fitMapqModelsLogistic <- function(x, incl.Xs=T, incl.repeats=F, incl.xd=F, incl.len=F, rescale=T, sca=1.2, maxit=60) {
	tab.all <- x
	opt <- bestSca(tab.all, maxit=maxit, rescale=rescale, sca=sca)
	sca <- opt$par
	fit.all <- modelAllLogistic(tab.all, incl.Xs=incl.Xs, incl.repeats=incl.repeats, incl.xd=incl.xd, incl.len=incl.len, replace.na=T, rescale=rescale, sca=sca)
	tab.all$model_mapq <- bt0and1(fit.all$mapq)
	tab.all$model_mapq2 <- bt0and1(fit.all$mapq2)
	return(list(tb=tab.all, fit=fit.all, err1=rankingError(tab.all, tab.all$model_mapq), err2=rankingError(tab.all, tab.all$mapq), rescale=rescale, sca=sca))
}

# Given filenames for SAM files emitted by two tools, analyze 
fitMapqModelsLogisticFile <- function(x, incl.Xs=T, incl.repeats=F, incl.xd=F, incl.len=F, rescale=T, sca=1.2) {
	tab <- openTabulatedSam(x)
	tab.all <- tab[selectAligned(tab),]
	return(fitMapqModelsLogistic(tab.all))
}

if(F) {
	require("multicore")
	lens = c(75, 100, 125, 150, 175)
	dir <- "/Users/langmead/cvs/bowtie2/paper/experiments/mason_bwa/"
	fn.all <- paste(dir, "r0_ill_", lens, "_100k.bt2_s.sat", sep="")
	tabl <- mclapply(fn.all,
		function(x) {
			xtab <- openTabulatedSam(x); xtab.al <- xtab[selectAligned(xtab),]; return(xtab.al)
		}, mc.cores=6)
	names(tabl) <- lens
	tabl$all <- do.call("rbind", tabl)
	
	# Fit each read length separately
	fitl <- mclapply(tabl,
		function(x) {
			fitMapqModelsLogistic(x)
		}, mc.cores=6)
	# Compare the mapping qualities from the overall fit with the qualities
	# from the fit over just the 100-length data.
	errs <- lapply(lens, function(x) {
		tab <- fitl[[toString(x)]]$tb;
		tabc <- fitl[["all"]]$tb[nchar(fitl[["all"]]$tb$seq) == x,];
		errc <- rankingError(tabc, tabc$model_mapq)
		err <- rankingError(tab, tab$model_mapq)
		return(c(errc, err, errc - err))
	})
}
