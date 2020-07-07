##
## PLOTS TWO LINES ON PHANEROZOIC AXES.
##

phanerozoicPlotBox <- function(y_top, y_bot, x_old, x_yng, y_lab) {

	source("/Volumes/shared/svnSifr/analyses/sharedCode/PhanerozoicBinDefs.R")
	
	y_rscale <-0
	y_lscale <-0

	# RANGE OF Y-AXIS
	if (missing(y_top)) {
		if (!missing(dataLeft)) 
			y_lscale <- max(dataLeft, na.rm=TRUE)
		if (!missing(dataRight)) 
			y_rscale <- max(dataRight, na.rm=TRUE)
		if (y_rscale < y_lscale) y_scale <- y_lscale
		if (y_rscale > y_lscale) y_scale <- y_rscale
		y_top <- y_scale*1.1
	} else {
		y_scale <- y_top
	}
	if (missing(y_bot))
		y_bot <- 0
	y_range <- y_top - y_bot
	
	if (missing(y_lab))
		y_lab <- ''
	
	# RANGE OF X-AXIS
	if (missing(x_old))
		x_old <- max(big)
	if (missing(x_yng))
		x_yng <- min(big)
	
	# TIME AXIS IS SET TO A PROPORTION OF Y-AXIS RANGE
	axis_height <- .1
	
	# LOCATION OF THE TOP OF THE TIME AXIS
	axis_top <- y_bot
	
	# SIZE OF LABLE TEXT
	time_cex <- 0.6
	
	h1 <- 0
	h2 <- (-1*axis_height*y_range/4)
	h3 <- (-1*axis_height*y_range/2) + h2
	h4 <- (-1*axis_height*y_range/2) + h3
			
	h1 <- h1 + axis_top
	h2 <- h2 + axis_top
	h3 <- h3 + axis_top
	h4 <- h4 + axis_top

	plot(1:1,type="n", ylim=c(h4,y_top), ylab=y_lab, xlim=c(x_old,x_yng), xlab="", main='', xaxs="i", xaxt="n", yaxs="i", yaxt="n",cex.axis=1.6, cex.main=1.6, cex.lab=1.6)

	# PREVENTS SEEING LINES PLOTTED BELOW AXIS
	polygon(c(x_old,x_yng,x_yng,x_old),c(h1,h1,h4,h4), col='white', bg='white')
		
	lines (c(x_old,x_yng),c(h1,h1))
	lines (c(x_old,x_yng),c(h2,h2))
	lines (c(x_old,x_yng),c(h3,h3))
	lines (c(x_old,x_yng),c(h4,h4))
	
	for (i in 1:length(bins)) {
		lines (c(bins[i],bins[i]),c(h2,h1))
	}
	
	for (i in 1:length(peri)) {
		lines (c(peri[i],peri[i]),c(h3,h2))
		xc <- (((peri[i+1]-peri[i])/2)+peri[i])
		yc <- h3+((h2-h3)/2)
		text (xc,yc,peri_n[i], cex=time_cex)
	}
	
	for (i in 1:length(big)) {
		lines (c(big[i],big[i]),c(h4,h3))
		xc <- (((big[i+1]-big[i])/2)+big[i])
		yc <- h4+((h3-h4)/2)
		text (xc,yc,big_n[i], cex=time_cex)
	}
}

