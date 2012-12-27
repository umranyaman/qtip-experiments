source("../shared.R")
source("../models/logistic.R")

if(F) {
	require("multicore")
	lens = c(75, 100, 125, 150, 175)
	dir <- "/Users/langmead/cvs/bowtie2/paper/experiments/mason_bwa/"
	fn.all <- paste(dir, "r0_ill_", lens, "_100k.bt2_s.sat", sep="")
	tabl <- mclapply(fn.all,
			function(x) {
				xtab <- openTabulatedSam(x); xtab.al <- xtab[selectAligned(xtab),]; return(xtab)
			}, mc.cores=6)
	names(tabl) <- lens
}

# An attempt to figure out the best way to include length in the model.  By
# using including alignment scores and differences between them in the model,
# we capture read length to an extent, since longer reads have wider ranges of
# scores than shorter reads.
#
# $err1
# $err1[[1]]
# [1] 10360516 10292714    67802
# 
# $err1[[2]]
# [1] 6023027 5977758   45269
# 
# $err1[[3]]
# [1] 4099201 4081545   17656
# 
# $err1[[4]]
# [1] 3007227 3014557   -7330
# 
# $err1[[5]]
# [1] 2394542 2390405    4137
# 
# 
# $err2
# $err2[[1]]
# [1] 10476561 10292663   183898
# 
# $err2[[2]]
# [1] 6042912 5978330   64582
# 
# $err2[[3]]
# [1] 4099580 4080398   19182
# 
# $err2[[4]]
# [1] 3019592 3014849    4743
# 
# $err2[[5]]
# [1] 2420784 2390214   30570
# 
exploreModelingLengths <- function() {
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
	
	evalModel <- function(fitl) {
		return(lapply(lens, function(x) {
			tab <- fitl[[toString(x)]]$tb;
			tabc <- fitl[["all"]]$tb[nchar(fitl[["all"]]$tb$seq) == x,];
			errc <- rankingError(tabc, tabc$model_mapq)
			err <- rankingError(tab, tab$model_mapq)
			return(c(errc, err, errc - err))
		}))
	}
	
	#
	# Strategy 1: let the alignment scores capture the read
	#
	
	# Fit each read length separately
	fitl.1 <- mclapply(tabl,
		function(x) {
			fitMapqModelsLogistic(x)
		}, mc.cores=6)
	# Compare the mapping qualities from the overall fit with the qualities
	# from the fit over just the 100-length data.
	errs.1 <- evalModel(fitl.1)
	
	#
	# Strategy 2: standardize the scores and include read length in model
	#
	
	# Fit each read length separately
	fitl.2 <- mclapply(tabl,
		function(x) {
			fitMapqModelsLogistic(x, incl.len=T, rescale=T)
		}, mc.cores=6)
	# Compare the mapping qualities from the overall fit with the qualities
	# from the fit over just the 100-length data.
	errs.2 <- evalModel(fitl.2)
	
	return(list(err1=errs.1, err2=errs.2))
}
