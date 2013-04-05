library(bitops)
library(zoo)

# Open a tabulated sam file
openTabulatedSam <- function(fn) {
	return(read.table(fn, comment.char="", quote="", header=T, stringsAsFactors=F))
}

# Given data frame of alignments, return vector of booleans indicating which aligned
selectAligned <- function(x) {
	return(bitAnd(x$flag, 4) == 0)
}

# Given data frame of alignments, return vector of booleans indicating which are unpaired
selectUnpaired <- function(x) {
	return(x$YT.Z == "UU")
}

# Given data frame of alignments, return vector of booleans indicating which are concordant
selectConcordant <- function(x) {
	return(x$YT.Z == "CP")
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

openTrainingAndSam <- function(x) {
	samTsvFn <- paste0(x, ".sat")
	trainingUnpFn <- paste0(x, ".unp.training.tsv")
	ret <- list(sam=read.table(samTsvFn, comment.char="", quote="", header=T, stringsAsFactors=F))
	# TODO: Split alignments up into unpaired, paired-end non-concordant, and paired-end concordant
	ret$al <- ret$sam[selectAligned(ret$sam),]
	ret$alUnp <- ret$al[selectUnpaired(ret$al),]
	ret$alConc <- ret$al[selectConcordant(ret$al),]
	if(file.exists(trainingUnpFn)) {
		ret$trainUnp <- read.table(trainingUnpFn, comment.char="", quote="", header=T, stringsAsFactors=F)
	}
	trainingDisc1Fn <- paste0(x, ".m1disc.training.tsv")
	trainingDisc2Fn <- paste0(x, ".m2disc.training.tsv")
	if(file.exists(trainingDisc1Fn)) {
		ret$trainDisc1 <- read.table(trainingDisc1Fn, comment.char="", quote="", header=T, stringsAsFactors=F)
		ret$trainDisc2 <- read.table(trainingDisc2Fn, comment.char="", quote="", header=T, stringsAsFactors=F)
	}
	trainingConc1Fn <- paste0(x, ".m1conc.training.tsv")
	trainingConc2Fn <- paste0(x, ".m2conc.training.tsv")
	trainingPairFn <- paste0(x, ".pair.training.tsv")
	if(file.exists(trainingConc1Fn)) {
		ret$trainConc1 <- read.table(trainingConc1Fn, comment.char="", quote="", header=T, stringsAsFactors=F)
		ret$trainConc2 <- read.table(trainingConc2Fn, comment.char="", quote="", header=T, stringsAsFactors=F)
		ret$trainPair <- read.table(trainingPairFn, comment.char="", quote="", header=T, stringsAsFactors=F)
	}
	return(ret)
}

alignmentsToTuples <- function(x) {
	diff <- as.numeric(x$AS.i - x$XS.i)
	diff <- ifelse(is.na(diff), 10000, diff)
	return(data.frame(
		len=nchar(x$seq),
		best=as.numeric(x$AS.i),
		diff=diff))
}

# Training data comes in the form of a data frame with columns

trainingMatrices <- function(dat) {
	trainUnp <- dat$trainUnp
	
	# Split training data up into repetitive, non-repetitive
	trainUnpR <- trainUnp[trainUnp$diff < 10000,]
	trainUnpU <- trainUnp[trainUnp$diff == 10000,]
	
	# Given a training-data data frame, turn it in to two matrices.
	minAs <- min(trainUnp$best)
	maxAs <- max(trainUnp$best)
	minDiff <- 0
	maxDiff <- max(trainUnpR$diff)
	
	factorAs <- factor(seq(maxAs, minAs))
	factorDiff <- factor(c(seq(minDiff, maxDiff)))
	
	trainUnpRCor <- trainUnpR[trainUnpR$correct == 1,]
	trainUnpRCor$best <- factor(trainUnpRCor$best, factorAs)
	trainUnpRCor$diff <- factor(trainUnpRCor$diff, factorDiff)
	
	trainUnpR$best <- factor(trainUnpR$best, factorAs)
	trainUnpR$diff <- factor(trainUnpR$diff, factorDiff)
	
	repCor <- as.matrix(table(trainUnpRCor$best, trainUnpRCor$diff))
	repTot <- as.matrix(table(trainUnpR$best, trainUnpR$diff))
	
	uniCor <- as.matrix(table(trainUnpU$best[trainUnpU$correct == 1]))
	uniTot <- as.matrix(table(trainUnpU$best))
	
	return(list(repCor=repCor, repTot=repTot, uniCor=uniCor, uniTot=uniTot))
}

knnFit <- function(x) {
	library(class)
	test <- alignmentsToTuples(x$alUnp)
	train <- x$trainUnp
	cl <- factor(train$correct)
	train$correct <- NULL
	return(knn(train, test, cl, 1, 0, prob=TRUE, use.all = TRUE))
}

compareRepeatCorrect <- function(x) {
	# Note, need to count both alignments where XS.i is defined *and*
	# alignments where Xs.i is defined.
	realCmp <- table(is.na(x$al$XS.i) & is.na(x$al$Xs.i), x$al$ZC.i == 1, dnn=c("Unique?", "Correct?"))
	simCmp <- table(x$trainUnp$diff == 10000, x$trainUnp$correct == 1, dnn=c("Unique?", "Correct?"))
}

# Given an table of aligned reads, return a list of vectors to use as
# alignment score-related covariates
getCovars <- function(x, incl.Xs=F, incl.xd=F, replace.na=F, rescale=F, sca=1.0) {
	if(replace.na) {
		repl <- max(x$AS.i) - sca * (max(x$AS.i) - min(x$AS.i))
		xs <- ifelse(is.na(x$XS.i), repl, x$XS.i)
		if(incl.Xs) {
			xs <- pmax(xs, x$Xs.i, na.rm=T)
		}
	} else {
		xs <- x$XS.i
	}
	xs <- x$AS.i - xs # could have some NAs
	as <- x$AS.i
	if(rescale) {
		as <- rescale(as, x$YN.i, x$Yn.i)
		xs <- rescale(xs, 0, sca * (x$Yn.i - x$YN.i))
	}
	return(list(as=as, xs=xs))
}

# Use Rafa's loess code to fit a loess model to the training data.
# Return a model that we can then use to predict
fitLoess <- function(dat, span=2.0, truncateK=T, truncateK.scale=1.5, scaleDiff="none") {
	repCor = dat$repCor 
	repTot = dat$repTot
	p = repCor / repTot
	tt = repTot
	x = as.numeric(rownames(p))
	y = as.numeric(colnames(p))
	#y[y == 10000] <- y[length(y)-1]+1
	#plot(x, rowMeans(p, na.rm=TRUE))
	#plot(y, colMeans(p, na.rm=TRUE))
	
	##after data exploration noticed that
	###if y>25 then it is pretty much 1, so let's get rid of them
	K <- max(y)
	if(truncateK) {
		K <- y[which(smooth(na.approx(colMeans(p, na.rm=TRUE))) >= 1.0)[1]]
		K <- as.integer(round(K * truncateK.scale))
		if(K < y[length(y)]) {
			repCor[, y == K] <- rowSums(repCor[, y >= K])
			repCor <- repCor[, y <= K]
			repTot[, y == K] <- rowSums(repTot[, y >= K])
			repTot <- repTot[, y <= K]
		}
	}
	p = repCor / repTot
	tt = repTot
	x = as.numeric(rownames(p))
	y = as.numeric(colnames(p))
	
	tmp = expand.grid(x, y)
	X = tmp[,1]
	Y = tmp[,2]
	if(scaleDiff == "log") {
		Y = log2(Y + 0.5)
	} else if(scaleDiff == "sqrt") {
		Y = sqrt(Y)
	}
	fit = loess(as.vector(p) ~ X * Y, span=span, degree=1, weight=sqrt(as.vector(tt)))
	
	#xbig <- seq(min(x), max(x))
	#ybig <- seq(min(y), max(y))
	#tmp = expand.grid(xbig, ybig)
	#Xbig = tmp[,1]
	#Ybig = tmp[,2]
	#z = matrix(predict(fit, newdata=data.frame(X=Xbig, Y=Ybig)), length(xbig), length(ybig))
	#dimnames(z) = list(as.character(xbig), as.character(ybig)) 
	#z[z > 1] <- 1
	#z[z < 0] <- 0
	
	###z is the smoothed map. upi can use predict above to create map
	
	#image(-x, y, z, col=rev(brewer.pal(9, "Blues")))
	
	###look at countour of function of x stratified by y
	###look at countour of function of y stratified by x
	#cols = rgb(seq(.9, .1, len=ncol(z)), 0, seq(.1, .9, len=ncol(z)))
	#matplot(x, z, col=cols, lty=1, type="l")
	#matplot(y, t(z), col=cols, lty=1, type="l")
	
	#return(list(z=z, maxDiff=K, fit=fit))
	return(list(maxDiff=K, fit=fit))
}

# Take table of alignments and parallel vector of mapping qualities and
# return a data frame containing columns for MAPQs (sorted in
# descending order), cumulative number of correct alignments, and
# cumulative number of incorrect alignments.
roc_table <- function(x, mapq) {
	xord <- order(mapq)
	tab <- table(mapq[xord], x$ZC.i[xord])
	tab <- data.frame(mapq=as.numeric(rownames(tab)), cor=tab[,2], incor=tab[,1])
	tab <- tab[order(tab$mapq, decreasing=T),]
	tab$corCum <- cumsum(tab$cor)
	tab$incorCum <- cumsum(tab$incor)
	return(tab)
}

plot_pcor <- function (data, pcor) {
	xord <- order(pcor)
	X <- pcor[xord]
	Y <- data$alUnp$ZC.i[xord]
	
}

# Plot a roc-like curve (except using absolute numbers of true
# positives and false negatives) for the given alignments (x) and
# mapping-quality ranking (mapq).
roc_plot <- function(x, mapq, color="blue", first=T, xlim=NULL, ylim=NULL, main="MAPQ ROCs", xlab="Incorrect alignments", ylab="Correct alignments", type="l") {
	tab <- roc_table(x, mapq)
	if(first) {
		plot(cumsum(tab$incor), cumsum(tab$cor), col=color, type=type, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main)
	} else {
		lines(cumsum(tab$incor), cumsum(tab$cor), col=color, type=type)
	}
}

# Plot two roc-like curves (except using absolute numbers of true
# positives and false negatives) for the given alignments (x) and the
# given two mapping-quality rankings (mapq1 and mapq2).
roc_2plot <- function(x, mapq1, mapq2, colors=c("blue", "red"), sca=1, round_digits=2, mapq1_lab="MAPQ 1", mapq2_lab="MAPQ 2", main="MAPQ ROCs") {
	if(round_digits > 0) {
		mapq1 <- round(mapq1, digits=round_digits)
		mapq2 <- round(mapq2, digits=round_digits)
	}
	xroc <- roc_table(x, mapq1)
	xroc2 <- roc_table(x, mapq2)
	mx_x <- max(max(cumsum(xroc$incor)), max(cumsum(xroc2$incor)))
	mn_x <- min(min(cumsum(xroc$incor)), min(cumsum(xroc2$incor)))
	mx_y <- max(max(cumsum(xroc$cor)), max(cumsum(xroc2$cor)))
	mn_y <- min(min(cumsum(xroc$cor)), min(cumsum(xroc2$cor)))
	roc_plot(x, mapq1,  color=colors[0], first=T, xlim=c(mn_x, mn_x + (mx_x - mn_x) * sca), ylim=c(mn_y + (mx_y - mn_y) * (1-sca), mx_y), main=main)
	roc_plot(x, mapq2, color=colors[1], first=F, xlim=c(mn_x, mn_x + (mx_x - mn_x) * sca), ylim=c(mn_y + (mx_y - mn_y) * (1-sca), mx_y))
	legend("bottomright", c(mapq1_lab, mapq2_lab), col=colors, pch=c(1, 1), lty=c(1, 1))
}

# Plot three roc-like curves (except using absolute numbers of true
# positives and false negatives) for the given alignments (x) and the
# given two mapping-quality rankings (mapq1 and mapq2).
roc_3plot <- function(x, mapq1, mapq2, mapq3, colors=c("blue", "red", "purple"), sca=1, round_digits=2, mapq1_lab="MAPQ 1", mapq2_lab="MAPQ 2", mapq3_lab="MAPQ 3", main="MAPQ ROCs") {
	if(round_digits > 0) {
		mapq1 <- round(mapq1, digits=round_digits)
		mapq2 <- round(mapq2, digits=round_digits)
		mapq3 <- round(mapq3, digits=round_digits)
	}
	xroc <- roc_table(x, mapq1)
	xroc2 <- roc_table(x, mapq2)
	xroc3 <- roc_table(x, mapq3)
	mx_x <- max(max(cumsum(xroc$incor)), max(cumsum(xroc2$incor)), max(cumsum(xroc3$incor)))
	mn_x <- min(min(cumsum(xroc$incor)), min(cumsum(xroc2$incor)), min(cumsum(xroc3$incor)))
	mx_y <- max(max(cumsum(xroc$cor)), max(cumsum(xroc2$cor)), max(cumsum(xroc3$cor)))
	mn_y <- min(min(cumsum(xroc$cor)), min(cumsum(xroc2$cor)), min(cumsum(xroc3$cor)))
	roc_plot(x, mapq1, color=colors[1], first=T, xlim=c(mn_x, mn_x + (mx_x - mn_x) * sca), ylim=c(mn_y + (mx_y - mn_y) * (1-sca), mx_y), main=main)
	roc_plot(x, mapq2, color=colors[2], first=F, xlim=c(mn_x, mn_x + (mx_x - mn_x) * sca), ylim=c(mn_y + (mx_y - mn_y) * (1-sca), mx_y))
	roc_plot(x, mapq3, color=colors[3], first=F, xlim=c(mn_x, mn_x + (mx_x - mn_x) * sca), ylim=c(mn_y + (mx_y - mn_y) * (1-sca), mx_y))
	legend("bottomright", c(mapq1_lab, mapq2_lab, mapq3_lab), col=colors, pch=c(1, 1, 1), lty=c(1, 1, 1))
}

# Penalty is sum of ranks of all incorrectly aligned reads.  Rank is lower for
# reads with lower MAPQ.
rankingError <- function(x, mapq) {
	ordr <- order(mapq)
	return(sum(as.numeric(which(x$ZC.i[ordr] == 0))))
}

predictLoess <- function(data, fit, scaleDiff="none") {
	xal <- data$alUnp
	cov <- getCovars(xal, incl.Xs=T, replace.na=T)
	as <- cov$as
	xs <- pmin(cov$xs, fit$maxDiff)
	if(scaleDiff == "log") {
		xs <- log2(xs + 0.5)
	} else if(scaleDiff == "sqrt") {
		xs <- sqrt(xs)
	}
	pcor <- predict(fit$fit, newdata=data.frame(X=as, Y=xs))
	pcor <- pmin(1.0, pmax(0.0, pcor))
	return(data.frame(pcor=pcor, mapq = -10 * log10(1.0 - pcor), as=as, xs=xs))
}

loadResults <- function(len, sensitivity, span=2.0, truncateK=F, truncateK.scale=1.5, scaleDiff="log") {
	fn <- paste0("r0_ill_mason_", len, "_50k.bt2_", sensitivity)
	print("Opening training and SAM...")
	data <- openTrainingAndSam(fn)
	print("Building matrices...")
	mats <- trainingMatrices(data)
	print("Fitting...")
	fit <- fitLoess(mats, span=span, truncateK=truncateK, truncateK.scale=truncateK.scale, scaleDiff=scaleDiff)
	print("Predicting...")
	prediction <- predictLoess(data, fit, scaleDiff=scaleDiff)
	pcor <- prediction$pcor
	mapq <- prediction$mapq
	print("Looking at performance...")
	err <- rankingError(data$alUnp, mapq)
	err.knn <- rankingError(data$alUnp, data$alUnp$XQ.f)
	err.old <- rankingError(data$alUnp, data$alUnp$mapq)
	roc <- roc_table(data$alUnp, round(mapq))
	roc.knn <- roc_table(data$alUnp, round(data$alUnp$XQ.f))
	roc.old <- roc_table(data$alUnp, round(data$alUnp$mapq))
	return(list(data=data, mats=mats, fit=fit, prediction=prediction, pcor=pcor, mapq=mapq, err=err, err.knn=err.knn, err.old=err.old, roc=roc, roc.knn=roc.knn, roc.old=roc.old))
}

if(FALSE) {
	setwd("~/Documents/workspace/mapq/examples/mason/")
	
	plotRoc <- function(results) {
		roc_3plot(results$data$alUnp, results$mapq, round(results$data$alUnp$XQ.f), results$data$alUnp$mapq, sca=0.5, mapq1_lab="Loess", mapq2_lab="Distance-weighted KNN", mapq3_lab="Bowtie 2")
	}
	
	plotRoc <- function(results) {
		roc_3plot(results$data$alUnp, results$mapq, results$data$alUnp$XQ.f, results$data$alUnp$mapq, sca=0.5, mapq1_lab="Loess", mapq2_lab="Distance-weighted KNN", mapq3_lab="Bowtie 2")
	}

	results.50  <- loadResults(50, "s")
	results.100 <- loadResults(100, "s")
	results.200 <- loadResults(200, "s")
	results.300 <- loadResults(300, "s")
	
	# Investigating horizontal jump
	# roc_table(results.300$data$alUnp, round(results.300$mapq, digits=1))
	# results.300.new$prediction[abs(results.300.new$prediction$mapq - 22.3) < 0.05,]
}

# We made a 'z' matrix using the training data

# We could make a corresponding 'z' matrix using the test data and then see where they differ the most
