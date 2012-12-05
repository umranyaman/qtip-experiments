source("shared.R")

# Take table of alignments and parallel vector of mapping qualities and
# return a data frame containing columns for MAPQs (sorted in
# descending order), cumulative number of correct alignments, and
# cumulative number of incorrect alignments.
roc_table <- function(x, mapq) {
	xord <- order(mapq)
	tab <- table(mapq[xord], x$ZC.i[xord])
	tab <- data.frame(mapq=as.numeric(rownames(tab)), cor=tab[,2], incor=tab[,1])
	tab <- tab[order(tab$mapq, decreasing=T),]
	return(tab)
}

# Return a 2-tuple containing the number of points in the x and y roc
# tables that are on the pareto frontier.  This is most intepretable if
# we divide the mapq spectra into the same number of strata.  By
# default, we divide them into 50 strata. 
roc_table_compare <- function(x, mapq1, mapq2, strata=50) {
	mapq1 <- bt0and1(round(bt0and1(mapq1)*strata))
	mapq2 <- bt0and1(round(bt0and1(mapq2)*strata))
	tab1 <- roc_table(x, mapq1)
	tab2 <- roc_table(x, mapq2)
	npareto <- c(0, 0)
	auc <- c(0, 0)
	cor1 <- cumsum(tab1$cor)
	cor2 <- cumsum(tab2$cor)
	incor1 <- cumsum(tab1$incor)
	incor2 <- cumsum(tab2$incor)
	for(i in 1:length(cor1)) {
		frontier <- T
		for(j in 1:length(cor2)) {
			if(cor2[j] >= cor1[i] && incor2[j] <= incor1[i]) {
				frontier <- F
				break
			}
		}
		if(frontier) {
			npareto[1] <- npareto[1] + 1
		}
	}
	for(i in 1:length(cor2)) {
		frontier <- T
		for(j in 1:length(cor1)) {
			if(cor1[j] >= cor2[i] && incor1[j] <= incor2[i]) {
				frontier <- F
				break
			}
		}
		if(frontier) {
			npareto[2] <- npareto[2] + 1
		}
	}
	return(npareto)
}

# 
roc_plot <- function(x, mapq, color="blue", first=T, xlim=NULL, ylim=NULL) {
	tab <- roc_table(x, mapq)
	if(first) {
		plot(cumsum(tab$incor), cumsum(tab$cor), col="blue", type="o", xlim=xlim, ylim=ylim)
	} else {
		lines(cumsum(tab$incor), cumsum(tab$cor), col="red", type="o", xlim=xlim, ylim=ylim)
	}
}

# 
roc_2plot <- function(x, mapq, mapq2, colors=c("blue", "red"), sca=1, round_digits=2) {
	if(round_digits > 0) {
		mapq <- round(mapq, digits=round_digits)
		mapq2 <- round(mapq2, digits=round_digits)
	}
	xroc <- roc_table(x, mapq)
	xroc2 <- roc_table(x, mapq2)
	mx_x <- max(max(cumsum(xroc$incor)), max(cumsum(xroc2$incor)))
	mn_x <- min(min(cumsum(xroc$incor)), min(cumsum(xroc2$incor)))
	mx_y <- max(max(cumsum(xroc$cor)), max(cumsum(xroc2$cor)))
	mn_y <- min(min(cumsum(xroc$cor)), min(cumsum(xroc2$cor)))
	print(mn_x)
	print(mx_x)
	print(mn_y)
	print(mx_y)
	roc_plot(x, mapq,  color=colors[0], first=T, xlim=c(mn_x, mn_x + (mx_x - mn_x) * sca), ylim=c(mn_y + (mx_y - mn_y) * (1-sca), mx_y))
	roc_plot(x, mapq2, color=colors[1], first=F, xlim=c(mn_x, mn_x + (mx_x - mn_x) * sca), ylim=c(mn_y + (mx_y - mn_y) * (1-sca), mx_y))
}
