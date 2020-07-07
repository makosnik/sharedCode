# NOTES
# FOR NORMALLY DISTRIBUTED DATA:
# http://en.wikipedia.org/wiki/Standard_deviation
# 1-SIGMA = 68.26894921371
# 2-SIGMA =	95.44997361036
# 3-SIGMA = 99.73002039367
# 4-SIGMA = 99.99366575163
#
# Chebyshev's inequality
#  At least 50% of the values are within 1.41 standard deviations from the mean
#  At least 75% of the values are within 2 standard deviations from the mean.
#  At least (1 - 1/k^2) x 100% of the values are within k standard deviations from the mean.

# MAKES HISTOGRAM WITH CI'S INDICATED
#
###################################################
histCI <- function(hdata, ci, bks, colBar, colCI, xlab, ylab, title) {

	if (missing(bks))
		bks <- 20

	if (missing(colBar))
		colBar <- 'grey'

	if (missing(colCI))
		colCI <- 'black'

	if (missing(xlab))
		xlab <- ''

	if (missing(ylab))
		ylab <- ''

	if (missing(title))
		title <- ''


	hd <- hist(hdata,col=colBar, breaks=bks, main=title, xlab=xlab, ylab=ylab)

	if (!missing(ci)) {
	
		hq <- quantile(hdata,ci,na.rm=TRUE)

		if (length(ci) > 1) {
			modY <- rep(c(0.05,-0.15),length.out=length(ci))
		} else {
			modY <- rep(0,2)
		}
		
		if (hq[length(hq)] < 1) {
			sigFig <- 3
		} else if (hq[length(hq)] < 10) {
			sigFig <- 2
		} else if (hq[length(hq)] < 100) {
			sigFig <- 1
		} else {
			sigFig <- 0
		}
		
		for (c in 1:length(ci)) {
			lines(c(hq[c],hq[c]),c(0,max(hd$counts)*(.78+modY[c])), lwd=1, lty=2, col=colCI)
			text(hq[c],max(hd$counts)*(.9+modY[c]),ci[c], cex=0.8)
			text(hq[c],max(hd$counts)*(.82+modY[c]),paste('(',round(hq[c],sigFig),')',sep=''), cex=0.8)
		}
	}
	
	return(hd)
}


###################################################
histCurve <- function(histData,xmin,xmax,ymax,bks,title,xlabel,ylabel,showBinSize,showRate,showRateCI,plotRateCI,REPS) {

	require(MASS)

	# SET - YEARS SINCE 1950
	yearCurrent <- 56

	# TRAP MISSING VALUES
	if (missing(xmin))
		xmin <- min(histData)

	if (missing(xmax))
		xmax <- max(histData)

	if (missing(title))
		title <- ''

	if (missing(xlabel))
		xlabel <- ''
	
	if (length(bks) == 1)
		bks <- seq(xmin,xmax+bks,bks)

	if (missing(ylabel)) {
		binSize <- round((bks[length(bks)]-bks[1])/length(bks),0)
		ylabel <- paste('shells per',binSize,'year bin')
	}
	
	if (missing(showBinSize))
		showBinSize <- FALSE

	if (missing(showRate))
		showRate <- TRUE

	if (missing(showRateCI))
		showRateCI <- TRUE

	if (missing(plotRateCI))
		plotRateCI <- TRUE

	if (missing(REPS))
		REPS <- 10000
	
	# NOT LIKE TO FIT EXPONENTIAL TO NEGATIVE X VALUES...
	MIN <- min(histData)
	if (MIN >= (-1*yearCurrent)) {
		MIN <- yearCurrent
	} else {
		MIN <- abs(MIN)
	}
	# HACK
	MIN <- 0
	
	# MAKE - HISTOGRAM OF AGES
	histogr <- hist(histData, breaks=bks, plot=F)

	# SET - YMAX
	if (missing(ymax)) {
		ymax <- max(histogr$density)
	} else {
		if (ymax < max(histogr$density)) ymax <- max(histogr$density)
	}
	
	# PLOT - HISTOGRAM
#	plot(histogr, freq=T, col='grey', xlab=xlabel, ylim=c(0,ymax), main=title, ylab=ylabel)
#	plot(histogr, freq=F, col='grey', xlab=xlabel, main=title, ylab=ylabel)

#	plotGeometric(histData, MIN, histogr, REPS, plotRateCI, showRate, showRateCI,title,xlabel,ylabel,xmax)
	plotExponential(histData, MIN, histogr, REPS, plotRateCI, showRate, showRateCI, title, xlabel, ylabel, xmax, ymax)
	
	if (showBinSize)
		text(xmax, ymax*0.7, pos=2, paste('bin size',bks,'yr'))

}

###################################################
plotGeometric <- function(histData, MIN, histogr, REPS, plotRateCI, showRate, showRateCI,title,xlabel,ylabel,xmax) {

	# SET - YMAX
	ymax <- max(histogr$density)

	#SET - GeoMids
	binMids <- round(histogr$mids+MIN,0)

	# FIT / MAKE - GEOMETRIC DISTRIBUTION TO AGES
	fitGeo <- fitdistr(histData+MIN, densfun='geometric')
	geoBest <- dgeom(binMids, prob=fitGeo$estimate)

	# CALCULATE CONFIDENCE INTERVALS - IF REQUIRED
	if (plotRateCI | showRateCI) {

		bootFit <- boot(histData+MIN, bootFitGeometric, R=REPS)
		bootCI <- boot.ci(bootFit, type="perc")

	}

	# MAKE / PLOT LINES BASED ON CONFIDENCE INTERVAL PARAMETERS
	if (plotRateCI) {

		geoBot <- dgeom(binMids, prob=bootCI$percent[4])
		geoTop <- dgeom(binMids, prob=bootCI$percent[5])

	}

	if (ymax < max(geoBest)) ymax <- max(geoBest)
	if (ymax < max(geoBot)) ymax <- max(geoBot)
	if (ymax < max(geoTop)) ymax <- max(geoTop)

	plot(histogr, freq=F, col='grey', ylim=c(0,ymax), xlab=xlabel, main=title, ylab=ylabel)
	lines(histogr$mids,geoBest, lty=1, lwd=2)
	lines(histogr$mids,geoBot, lty=2, lwd=1)
	lines(histogr$mids,geoTop, lty=2, lwd=1)


	# PLOT TEXT (x 10^-3 to make it easier to compare)
	plotText <- 1.0;
	if (showRate)
		text(xmax, ymax*0.9, pos=2, cex=plotText, paste('prob =',round(fitGeo$estimate*1000,2),'kyr'))
	if (showRateCI)
		text(xmax, ymax*0.8, pos=2, cex=plotText, paste('95% =',round(bootCI$percent[4]*1000,1),'-',round(bootCI$percent[5]*1000,1),'kyr'))

}

###################################################
plotExponential <- function(histData, MIN, histogr, REPS, plotRateCI, showRate, showRateCI, title, xlabel, ylabel, xmax, ymax) {	

	require(boot)
	source('../sharedCode/bootstrap.R')

	if (missing(xmax))
		xmax <- max(histogr$breaks)
	
	# SET - YMAX
	if (missing(ymax)) {
		# HACK FOR 1st time averging paper	
		# ymax <- 0.0031

		# HACK FOR 1st time averging paper	
		# if (ymax > 0.01) ymax <- 0.019
		# if (ymax < max(histogr$density)) ymax <- max(histogr$density)

		ymax <- 0.0001
	}
	
	if (ymax < max(histogr$density)) ymax <- max(histogr$density)

	#SET - Bin mid points
	binMids <- round(histogr$mids+MIN,0)

	# FIT / MAKE - DISTRIBUTION TO AGES
	fitDist <- fitdistr(histData+MIN, densfun='exponential')
	disBest <- dexp(binMids, rate=fitDist$estimate)

	# CALCULATE CONFIDENCE INTERVALS - IF REQUIRED
	if (plotRateCI | showRateCI) {

		bootFit <- boot(histData+MIN, bootFitExponential, R=REPS)
		bootCI <- boot.ci(bootFit, type="perc")

	}

	# MAKE / PLOT LINES BASED ON CONFIDENCE INTERVAL PARAMETERS
	if (plotRateCI) {

		disBot <- dexp(binMids, rate=bootCI$percent[4])
		disTop <- dexp(binMids, rate=bootCI$percent[5])

	}

	# SET Y-MAX
	if (ymax < max(disBest)) ymax <- max(disBest)
	if (ymax < max(disBot)) ymax <- max(disBot)
	if (ymax < max(disTop)) ymax <- max(disTop)

#	print(paste(title,'pe -',ymax))


	cMax <- max(histogr$counts)
	dMax <- max(histogr$density)
	vOne <- dMax/cMax

	if (cMax == 13) cMax <- 12
	if (cMax == 21) cMax <- 20
	if (cMax == 9) cMax<- 8
	if (cMax == 5) cMax<- 8
	if (cMax < 3) cMax <- 3

	aLen <- 5
	if (aLen > cMax) aLen <- cMax+1

	aNam <- seq(0,cMax,length.out=aLen)
	aLoc <- seq(0,(cMax*vOne),length.out=aLen)

##	13
##	22

	# PLOT HISTOGRAM WITH LINES
#	plot(histogr, freq=F, col='grey', xlab=xlabel, main=title, ylab=ylabel)
	plot(histogr, freq=F, col='grey', ylim=c(0,ymax), xlab=xlabel, main=title, ylab=ylabel,yaxt='n')
#	plot(histogr, freq=F, col='grey', ylim=c(0,ymax), xlab=xlabel, main=title, ylab=ylabel)
	axis(2, at=aLoc, labels=aNam, par(las=1))
	lines(histogr$mids,disBest, lty=1, lwd=2)
	lines(histogr$mids,disBot, lty=2, lwd=1)
	lines(histogr$mids,disTop, lty=2, lwd=1)

	# CALC
	
	#
	# mean life = 1 / fitExp$estimate
	# half life = ln(2) /  fitExp$estimate
	
	value <- 'half life'
	vBest <- round(log(2)/fitDist$estimate,0) 
	vMed <- round(log(2)/quantile(bootFit$t,probs=c(0.5)),0) 
	vBot <- round(log(2)/bootCI$percent[5],0) 
	vTop <- round(log(2)/bootCI$percent[4],0) 


	# PLOT TEXT (x 10^-3 to make it easier to compare)
	plotText <- 1.0;
	if (showRate)
		text(xmax, ymax*0.9, pos=2, cex=plotText, paste(value,'=',vBest,'yr'))
	if (showRateCI)
		text(xmax, ymax*0.8, pos=2, cex=plotText, paste('95% =',vBot,'-',vTop,'yr'))

#	text(xmax, ymax*0.6, pos=2, cex=plotText, paste('y/c/dmax =',round(ymax,4),'-',round(cMax,4),'-',round(dMax,4),'-',aLen))

}

###################################################
histCurveJack <- function(histData,xmin,xmax,bks,title,xlabel,ylabel) {

	# SET - YEARS SINCE 1950
	yearCurrent <- 56

	# TRAP MISSING VALUES
	if (missing(xmin))
		xmin <- min(histData)

	if (missing(xmax))
		xmax <- max(histData)

	if (missing(title))
		title <- ''

	if (missing(xlabel))
		xlabel <- ''

	if (missing(ylabel))
		ylabel <- ''

	# NOT LIKE TO FIT EXPONENTIAL TO NEGATIVE X VALUES...
	MIN <- min(histData)
	if (MIN >= (-1*yearCurrent)) {
		MIN <- yearCurrent
	} else {
		MIN <- abs(MIN)
	}
	
	# MAKE - HISTOGRAM OF AGES
	histogr <- hist(histData, breaks=seq(xmin,xmax+bks,bks), plot=F)

	# SET - YMAX
	ymax <- max(histogr$counts)
		
	# PLOT - HISTOGRAM
	plot(histogr, freq=T, col='grey', ylab=ylabel, xlab=xlabel, ylim=c(0,ymax), main=title)

	# SETUP - ARRAY TO STORE RESAMPLED FITS
	resArra <- array(data=0,length(histData))

	# RESAMPLE FIT - EXPONENTIAL DISTRIBUTION TO AGES
	for (i in 1:length(histData)) {
		selArra <- array(data=TRUE,length(histData))
		selArra[i] <- FALSE
		fitData <- subset(histData, selArra)
		resArra[i] <- fitdistr(fitData+MIN, densfun='exponential')$estimate
	}
	
	# SET - CONFIDENCE INTERVALS
	ci <- quantile(resArra,c(0.975,0.5,0.025))

	# CALCULATE CI DISTRIBUTIONS
	for (rate in ci) {
	
		# MAKE - EXPONENTIAL DISTRIBUTION TO PLOT
		x <- dexp(histogr$mids+MIN, rate=rate)

		# PLOT - FIT EXPONENTIAL
		lines(histogr$mids,x*max(histogr$counts)*1000, lty=3, lwd=3)
	}
	

	# PLOT TEXT (INSTANTANEOUS RATE)
	ci <- ci*1000
	text((xmax), (ymax*0.1)*9, pos=2, paste('median rate =',round(ci[2],2),'x',expression(10^-3)))
	text((xmax), (ymax*0.1)*8, pos=2, paste('(',round(ci[3],2),' - ',round(ci[1],2),')',sep=''))
	text((xmax), (ymax*0.1)*7, pos=2, paste('bin size',bks,'yr'))

}

