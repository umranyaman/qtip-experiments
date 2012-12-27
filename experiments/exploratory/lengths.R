source(file.path('..', 'shared.R'), chdir=T)
source(file.path('..', 'models', 'logistic.R'), chdir=T)

# An attempt to figure out the best way to include length in the model.  By
# using including alignment scores and differences between them in the model,
# we capture read length to an extent, since longer reads have wider ranges of
# scores than shorter reads.
# 
# > errs.1
# [[1]]
# [1] 10339087 10281300    57787
# 
# [[2]]
# [1] 6006029 5976588   29441
# 
# [[3]]
# [1] 4099877 4079940   19937
# 
# [[4]]
# [1] 3004331 3010722   -6391
# 
# [[5]]
# [1] 2385858 2351157   34701
# 
# > errs.2
# [[1]]
# [1] 10358917 10281300    77617
# 
# [[2]]
# [1] 5985028 5976588    8440
# 
# [[3]]
# [1] 4086075 4079940    6135
# 
# [[4]]
# [1] 3032209 3010722   21487
# 
# [[5]]
# [1] 2348697 2351157   -2460
# 
exploreModelingLengths <- function() {
	lens = c(75, 100, 125, 150, 175)
	dir <- "/Users/langmead/cvs/bowtie2/paper/experiments/mason_bwa/"
	fn.all <- paste(dir, "r0_ill_", lens, "_100k.bt2_s.sat", sep="")
	tabl <- lapply(fn.all,
			function(x) {
				xtab <- openTabulatedSam(x); xtab.al <- xtab[selectAligned(xtab),]; return(xtab.al)
			})
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
	
	# Print warnings as they occur
	options(warn=1)
	
	#
	# Strategy 1: let the alignment scores capture the read
	#
	
	# Fit each read length separately
	fitl.1 <- lapply(tabl,
		function(x) {
			fitMapqModelsLogistic(x)
		})
	# Compare the mapping qualities from the overall fit with the qualities
	# from the fit over just the 100-length data.
	errs.1 <- evalModel(fitl.1)
	
	#
	# Strategy 2: standardize the scores and include read length in model
	#
	
	# Fit each read length separately
	fitl.2 <- lapply(tabl,
		function(x) {
			fitMapqModelsLogistic(x, incl.len=T, rescale=T)
		})
	# Compare the mapping qualities from the overall fit with the qualities
	# from the fit over just the 100-length data.
	errs.2 <- evalModel(fitl.2)
	
	save(fitl.1, errs.1, fitl.2, errs.2, file="lengths.rda")
}

exploreModelingLengths()
