source("shared.R")

# Given a data frame containing all aligned reads (including those that
# align just once or many times) return a fit for a linear model with
# response = mapq ranking and terms = alignment scores.
#
# Second argument is a scaling factor, which determines how XS.i values
# should be assigned to reads that align only once.  Right now we set
# them equal to max(AS.i) - 2 * (max(AS.i) - min(AS.i)).  
modelLoess <- function(x, incl.Xs=F, incl.repeats=F, incl.xd=F, replace.na=T, sca=2.0) {
	cov <- getCovars(x, incl.Xs=incl.Xs, incl.xd=incl.xd, replace.na=replace.na, sca=sca)
	form <- "x$ZC.i ~ cov$as * cov$xs"
	if(incl.repeats) { form <- paste(form, "+ x$AllRepeats") }
	if(incl.xd) { form <- paste(form, "+ x$XD.i") }
	form <- as.formula(form)
	model <- loess(form)
	mapq <- predict(model, data.frame(as=cov$as, xs=cov$xs))
	return(list(model=model, mapq=mapq, coef=coef))
}

# Find the best scaling factor
bestScaLoess <- function(tab, incl.Xs=F, incl.repeats=F, incl.xd=F, sca=1.2, maxit=25, method="SANN") {
	fn <- function(x) {
		fit <- modelLoess(
				tab,
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
fitMapqModelsLoess <- function(x, incl.Xs=T, incl.repeats=F, incl.xd=F, sca=1.2) {
	tab.all <- x
	opt <- bestScaLoess(tab.all, method="Nelder-Mead", maxit=1, sca=sca)
	sca <- opt$par
	fit.all <- modelLoess(tab.all, incl.Xs=incl.Xs, incl.repeats=incl.repeats, incl.xd=incl.xd, replace.na=T, sca=sca)
	tab.all$model_mapq <- bt0and1(fit.all$mapq)
	return(list(tb=tab.all, fit=fit.all, err1=rankingError(tab.all, tab.all$model_mapq), err2=rankingError(tab.all, tab.all$mapq)))
}

# Given filenames for SAM files emitted by two tools, analyze 
fitMapqModelsLoessFile <- function(x, incl.Xs=T, incl.repeats=F, incl.xd=F, sca=1.2) {
	tab <- openTabulatedSam(x)
	tab.all <- tab[selectAligned(tab),]
	return(fitMapqModelsLoess(tab.all, incl.Xs=incl.Xs, incl.repeats=incl.repeats, incl.xd=incl.xd, sca=sca))
}

if(F) {
	dr <- "/Users/langmead/Documents/workspace/mapq"
	x <- paste(dr, "examples/mason/r0_ill_100_100k.bt2_s.sat.bz2", sep="/")
	fit <- fitMapqModelsLoessFile(x)
	roc_table_compare(fit$tb, fit$tb$model_mapq, fit$tb$mapq)
}
