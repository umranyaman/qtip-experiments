source("shared.R")
source("logistic.R")

# Plot three roc-like curves, one for a competitor of Bowtie 2, and two
# for Bowtie 2 - one using model mapping qualities and the other using
# native mapping qualities. 
cmp_plot <- function(x, xmapq1, xmapq2, y, ymapq, colors=c("blue", "red", "darkpink4"), sca=1, round_digits=2) {
	if(round_digits > 0) {
		xmapq1 <- round(xmapq1, digits=round_digits)
		xmapq2 <- round(xmapq2, digits=round_digits)
	}
	xroc1 <- roc_table(x, xmapq1)
	xroc2 <- roc_table(x, xmapq2)
	yroc <- roc_table(y, ymapq)
	mx_x <- max(max(cumsum(xroc$incor)), max(cumsum(xroc2$incor)))
	mn_x <- min(min(cumsum(xroc$incor)), min(cumsum(xroc2$incor)))
	mx_y <- max(max(cumsum(xroc$cor)), max(cumsum(xroc2$cor)))
	mn_y <- min(min(cumsum(xroc$cor)), min(cumsum(xroc2$cor)))
	roc_plot(x, mapq1,  color=colors[0], first=T, xlim=c(mn_x, mn_x + (mx_x - mn_x) * sca), ylim=c(mn_y + (mx_y - mn_y) * (1-sca), mx_y))
	roc_plot(x, mapq2, color=colors[1], first=F, xlim=c(mn_x, mn_x + (mx_x - mn_x) * sca), ylim=c(mn_y + (mx_y - mn_y) * (1-sca), mx_y))
	legend("bottomright", c("MAPQ 1", "MAPQ 2"), col=colors, pch=c(1, 1), lty=c(1, 1))
}

if(F) {
	dr <- "/Users/langmead/Documents/workspace/mapq"
	bt2 <- paste(dr, "examples/mason/r0_ill_100_100k.bt2_s.sat.bz2", sep="/")
	bwa <- paste(dr, "examples/mason/r0_ill_100_100k.bwa.sat.bz2", sep="/")
	bt2.tab <- openTabulatedSam(bt2)
	bwa.tab <- openTabulatedSam(bwa)
	bt2.al <- bt2.tab [selectAligned(bt2.tab ),]
	bwa.al <- bwa.tab [selectAligned(bwa.tab ),]
	bt2.fit <- fitMapqModelsLogistic(bt2.al)
	cmp_plot(bt2.al, fit$tb$model_mapq, bt2.al$mapq)
}
