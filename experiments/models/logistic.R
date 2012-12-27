source(file.path('..', 'shared.R'), chdir=T)

# Given a data frame containing all aligned reads (including those that
# align just once or many times) return a fit for a logistic regression
# model with response = mapq ranking and terms = alignment scores.
#
# "sca" is a scaling factor, which determines how XS.i values should be
# assigned to reads that align only once.  Right now we set them equal
# to max(AS.i) - 2 * (max(AS.i) - min(AS.i)).  
modelAllLogistic <- function(
	x,              # alignment table - only aligned reads
	incl.Xs=T,      # use Xs.i field when calculating score gap
	incl.repeats=F, # include AllRepeats as term in model
	incl.xd=F,      # include XD.i as term in model
	incl.len=F,     # include read length as term in model
	replace.na=T,   # replace NA XS.i field with gaps = (sca * max gap)
	rescale=F,      # standardize alignment scores, differences to be in [0, 1]
	sca=1.0)        # scaling factor to use when replacing NA XS.i field
{
	# Get vectors for main covariates (scores and gaps)
	cov <- getCovars(
		x,
		incl.Xs=incl.Xs,
		incl.xd=incl.xd,
		replace.na=replace.na,
		rescale=rescale,
		sca=sca)
	if(any(is.na(cov$as))) {
		warning("Error: one or more AS:is are NA")
	}
	if(length(cov$as) != length(cov$xs)) {
		warning("Error: best score vector and gap vector are different lengths")
	}
	form <- "x$ZC.i ~ cov$as + cov$xs"
	# Add extra covariates if requested
	if(incl.repeats && length(unique(nchar(x$AllRepeats))) > 1) {
		form <- paste(form, "+ x$AllRepeats")
	}
	if(incl.xd && length(unique(nchar(x$XD.i))) > 1) {
		form <- paste(form, "+ x$XD.i")
	}
	if(incl.len && length(unique(nchar(x$seq))) > 1) {
		form <- paste(form, "+ nchar(x$seq)")
	}
	# Fit logistic regression model
	model <- glm(as.formula(form), family=binomial("logit"))
	# Extract coefficients
	coef <- coefficients(model)
	# Calculate predicted p's
	ex <- coef[1] + (coef[2] * cov$as) + (coef[3] * cov$xs)
	next_coef <- 4
	if(incl.repeats && length(unique(nchar(x$AllRepeats))) > 1) {
		ex <- ex + coef[next_coef] * x$AllRepeats; next_coef <- next_coef + 1
	}
	if(incl.xd && length(unique(nchar(x$XD.i))) > 1) {
		ex <- ex + coef[next_coef] * x$XD.i; next_coef <- next_coef + 1
	}
	if(incl.len && length(unique(nchar(x$seq))) > 1) {
		ex <- ex + coef[next_coef] * nchar(x$seq); next_coef <- next_coef + 1
	}
	mapq <- (1.0 / (1.0 + exp(-ex)))
	if(any(is.na(mapq))) {
		warning("Error: one or more MAPQs are NA")
		print(summary(model))
		print(coef)
	}
	return(list(model=model, mapq=mapq, coef=coef))
}

# Use "optim" to search for the best scaling factor to use for the dataset
bestSca <- function(
	tab,            # alignment table - only aligned reads
	incl.Xs=T,      # use Xs.i field when calculating score gap
	incl.repeats=F, # include AllRepeats as term in model
	incl.xd=F,      # include XD.i as term in model
	incl.len=F,     # include read length as term in model
	rescale=F,      # standardize alignment scores, differences to be in [0, 1]
	sca=1.0)        # how to substitute 
{
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
		print(paste("Called with x=", x, " err=", err))
		return(err)
	}
	opt <- optimize(fn, c(0.1, 1.5), tol = 0.001)
	return(opt$minimum)
}

# Given filenames for SAM files emitted by two tools, analyze 
fitMapqModelsLogistic <- function(
	x,
	incl.Xs=T,
	incl.repeats=F,
	incl.xd=F,
	incl.len=F,
	rescale=F,
	sca=1.0)
{
	tab.all <- x
	sca <- bestSca(tab.all, rescale=rescale, sca=sca)
	fit.all <- modelAllLogistic(
		tab.all,
		incl.Xs=incl.Xs, incl.repeats=incl.repeats, incl.xd=incl.xd,
		incl.len=incl.len, replace.na=T, rescale=rescale, sca=sca)
	tab.all$model_mapq <- bt0and1(fit.all$mapq)
	return(list(tb=tab.all, fit=fit.all, err1=rankingError(tab.all, tab.all$model_mapq), err2=rankingError(tab.all, tab.all$mapq), rescale=rescale, sca=sca))
}

# Given filenames for SAM files emitted by two tools, analyze 
fitMapqModelsLogisticFile <- function(x, incl.Xs=T, incl.repeats=F, incl.xd=F, incl.len=F, rescale=F, sca=1.0) {
	tab <- openTabulatedSam(x)
	tab.all <- tab[selectAligned(tab),]
	return(fitMapqModelsLogistic(tab.all, incl.Xs=incl.Xs, incl.repeats=incl.repeats, incl.xd=incl.xd, incl.len=incl.len, rescale=rescale, sca=sca))
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
	fitl.1 <- mclapply(tabl,
		function(x) {
			fitMapqModelsLogistic(x)
		}, mc.cores=6)
	# Compare the mapping qualities from the overall fit with the qualities
	# from the fit over just the 100-length data.
	errs.1 <- lapply(lens, function(x) {
		tab <- fitl[[toString(x)]]$tb;
		tabc <- fitl[["all"]]$tb[nchar(fitl[["all"]]$tb$seq) == x,];
		errc <- rankingError(tabc, tabc$model_mapq)
		err <- rankingError(tab, tab$model_mapq)
		return(c(errc, err, errc - err))
	})

	# Fit each read length separately
	fitl.2 <- mclapply(tabl,
		function(x) {
			fitMapqModelsLogistic(x)
		}, mc.cores=6)
}
