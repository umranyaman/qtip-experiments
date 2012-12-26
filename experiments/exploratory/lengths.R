source("../shared.R")

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
