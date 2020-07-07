identifyOutliers <- function(qData, ci, qLabel, sLabel, xLabel, yLabel) {

	# GET QUANTILES
	hq <- quantile(qData[,qLabel],ci)

	# FOR EACH QUANTILE
	for (c in length(ci):1) {
		
		# GET LIST OF OUTLIERS
		tmp <-sumfd[qData[,qLabel] > hq[c],]
		if (is.matrix(bob)) {
			specN <- min(tmp[,sLabel])
		} else {
			specN <- tmp[sLabel]
		}
		
		specs <- rMatrix[rMatrix[,'n']>=specN,'specimen_no']
		
		for (s in 1:length(specs)) {
			xVal <- datesTaxon[(datesTaxon[,'specimen_no']==specs[s]),colAges[1]]^expTax
			yVal <- datesTaxon[(datesTaxon[,'specimen_no']==specs[s]),colAges[2]]
			points(xVal,yVal,pch=21, bg=ciCol[c], col=ciCol[c], cex=0.75, lwd=0.5)
			if (c==1) text(xVal,yVal,specs[s],cex=0.8,pos=1)
		}
	}
	
	return (hq)
}


#
#	REGRESSION PLOT - W/ ERROR BARS
# low varaince  ~ high variance
# predicted ~ predictor
###################################################
regressPlot <- function(x, y, xl, yl, ti, pr2, peb, pun, pul, plab, specNames, cCut) {

	if (missing(xl))
		xl<- ''

	if (missing(yl))
		yl<- ''

	if (missing(ti))
		ti<- ''

	if (missing(pr2))
		pr2 <- FALSE

	if (missing(peb))
		peb <- FALSE
	if (!is.matrix(y))
		peb <- FALSE
	if (!is.matrix(x))
		peb <- FALSE

	if (missing(pun))
		pun <- FALSE

	if (missing(pul))
		pul <- FALSE

	if (missing(plab))
		plab <- FALSE

	if (missing(cCut))
		cCut <- 0.99

	xMin <- min(x, na.rm=TRUE)
	xMax <- max(x, na.rm=TRUE)
	yMin <- min(y, na.rm=TRUE)
	yMax <- max(y, na.rm=TRUE)
	
	# IF DOING UNITY LINES, MAKE SURE THE PLOT IS SQUARE
	if (pul) {
		if (xMax > yMax) mMax <- xMax else mMax <- yMax
		if (xMin < yMin) mMin <- xMin else mMin <- yMin
		
		xMax <- yMax <- mMax
		xMin <- yMin <- mMin
		
	 	plot(1:1, type='n', ylab=yl, xlab=xl, main=ti, xlim=c(mMin,mMax), ylim=c(mMin,mMax))
	} else {
	 	plot(1:1, type='n', ylab=yl, xlab=xl, main=ti, xlim=c(xMin,xMax), ylim=c(yMin,yMax))
	}
	
	# PLOT POINTS
	
	if (is.matrix(y)) {
	 	r <- lm(y[,'fit'] ~ x[,'fit'] + 0)
		points(x[,'fit'], y[,'fit'], col='black', pch=22, cex=1)
		n <- length(x[,'fit'])
	} else {
	 	r <- lm(y ~ x)
		points(x, y, pch=21, col='white', bg='black', cex=1)
		n <- length(x)
	}
 	abline(r)

	# LABEL SPECIMENS WITH RESIDUALS OUTSIDE CONFIDENCE INTERVAL
	if ((plab) & (!missing(specNames))) {

		pData <- cbind(x, y, abs(r$residuals))

		# GET OFFSET FOR TEXT
		xmin <- min(pData[,'y'], na.rm=TRUE)
		xmax <- max(pData[,'y'], na.rm=TRUE)
		xoff <- (xmax - xmin)*0.05
	
		hq <- quantile(pData[,3],cCut)

		# PLOT SPECIMEN NUMBERS FOR THE REALLY BAD SPECIMENS 
		for (p in 1:length(pData[,'y'])) {
			if (pData[p,3]>hq[1]) {
				text(pData[p,'x'],pData[p,'y']-xoff,specNames[p], cex=0.8)
				points(pData[p,'x'],pData[p,'y'], pch=21, col='grey', bg='grey', cex=0.75, lwd=0.25)
			}
		}
	}

	# PLOT R2 VALUE
 	if (pr2) {
 		fontSize <- 1.0
		text(xMin,(yMin+(yMax-yMin)*0.9),cex=fontSize,pos=4,substitute(r^2 == k, list(k = round(summary(r)$adj.r.squared,3))))
		text(xMin,(yMin+(yMax-yMin)*0.8),cex=fontSize,pos=4,substitute(n == k, list(k = n)))
	}
	
	# PLOT ERROR BARS
	if (peb) {
		for (i in 1:length(x[,'fit'])) {
			lines(c(x[i,'fit'],x[i,'fit']),c(y[i,'upr'],y[i,'lwr']), lwd=1.0)
			lines(c(x[i,'upr'],x[i,'lwr']),c(y[i,'fit'],y[i,'fit']), lwd=1.0)
		}
		points(x[,'fit'], y[,'fit'], col='black', pch=22, cex=0.8, lwd=1.0)
	}
	# PLOT UNITY LINE
	if (pul) {
		lines (c(mMin,mMax),c(mMin,mMax), lwd=2, lty=2)
	}

	
	# PLOT UNCERTAINTY
	if (pun) {
		r.pred <- predict(r, int='p')
		r.var <- median(r.pred[,'upr'] - r.pred[,'lwr'])
		text(xMin,(yMin+(yMax-yMin)*0.7),cex=fontSize,pos=4,substitute(unc == k, list (k = round(r.var,0))))

		# PLOT - ADD CONFIDENCE AND PREDICTION INTERVALS
#		matlines(x,r.pred, lty=c(1,1,1), col=c('black','grey66', 'grey66'), lwd=2)
	}

	# RETURN PREDICTION INTERVAL	
	if (pun)
		return (r.pred)
}


###################################################
optimExpX <- function(exps, ...) {
	
	load('tmp.Rimage')

	fit <- lm(y ~ I(x^exps))
	
	return(1-summary(fit)$adj.r.squared)
}


#
#	REGRESSION PLOT - W/ ERROR BARS
# low varaince  ~ high variance
# predicted ~ predictor
#
# x, y		 = x & y data
# xl, yl, ti = x & y axis labels & title
# pl		 = sample numbers to plot on points
# pr2 F/T	 = plot r2, fit exponent and number of specimens
# peb F/T	 = plot error bars
# pun F/T	 = plot uncertanty
# pul F/T	 = plot unity line
###################################################
regressPlotOpim <- function(x, y, xl, yl, ti, pl, pr2, peb, pun, pul, plab, specNames, cCut, xLabExp, upLegend) {

#	library(lmodel2)

	if (missing(xl))
		xl<- ''

	if (missing(yl))
		yl<- ''

	if (missing(ti))
		ti<- ''

	if (missing(pr2))
		pr2 <- FALSE

	if (missing(peb))
		peb <- FALSE

	if (missing(pun))
		pun <- FALSE

	if (missing(pul))
		pul <- FALSE

	if (missing(plab))
		plab <- FALSE

	if (missing(cCut))
		cCut <- 0.99
	
	if (missing(xLabExp))
		xLabExp <- FALSE

	if (missing(upLegend))
		upLegend <- FALSE

	save(x, y, file ='tmp.Rimage')

	# OPTIMIZE EXPONENT FIT TO DATA
	optPara <- optimize(optimExpX, interval = c(0,6))
	j <- optPara$minimum

	xMin <- min(x^j, na.rm=TRUE)
	xMax <- max(x^j, na.rm=TRUE)
	yMin <- min(y, na.rm=TRUE)
	yMax <- max(y, na.rm=TRUE)
	
	# IF DOING UNITY LINES, MAKE SURE THE PLOT IS SQUARE
	if (pul) {
		if (xMax > yMax) mMax <- xMax else mMax <- yMax
		if (xMin < yMin) mMin <- xMin else mMin <- yMin
	}
	
	alab <- round(seq(min(x),max(x), length.out=5),2)
	aloc <- alab^j
	
	
	# PLOT POINTS
 	plot(1:1, type='n', ylab=yl, xlab='', main=ti, xlim=c(xMin,xMax), ylim=c(yMin,yMax), xaxt='n')
  
	if (xLabExp) {
		axis(side=1, at=aloc, labels=alab)
		title(xlab=xl)
	} else {
		axis(side=1, at=alab, labels=alab)
		xla <- substitute(x ^ k, list(x = xl, k = round(j,3)))
		title(xlab=xla)
 	}
 	
 	r <- lm(y ~ I(x^j))
 	abline(r)

	if(r$coefficients['I(x^j)'] > 0) {
		upLegend <- FALSE
	} else {
		upLegend <- TRUE
	}

	points(((x^j)),((y)), pch=21, col='white', bg='black', cex=0.75, lwd=0.25)

	# LABEL SPECIMENS WITH RESIDUALS OUTSIDE CONFIDENCE INTERVAL
	if ((plab) & (!missing(specNames))) {

		pData <- cbind(x, y, specNames, abs(r$residuals))

		# GET OFFSET FOR TEXT
		xmin <- min(pData[,'y'], na.rm=TRUE)
		xmax <- max(pData[,'y'], na.rm=TRUE)
		xoff <- (xmax - xmin)*0.05
	
		hq <- quantile(pData[,4],cCut)

		# PLOT SPECIMEN NUMBERS FOR THE REALLY BAD SPECIMENS 
		for (p in 1:length(pData[,'y'])) {
			if (pData[p,4]>hq[1]) {
				text(pData[p,'x']^j,pData[p,'y']-xoff,specNames[p], cex=0.8)
				points(pData[p,'x']^j,pData[p,'y'], pch=21, col='grey', bg='grey', cex=0.75, lwd=0.25)
			}
		}
	}

	# PLOT R2 VALUE
 	if (pr2) {
 		fontSize <- 1.0
  		if (upLegend) {
			xPo <- min(x^j)
		} else {
			xPo <- (min(x^j) + (max(x^j) - min(x^j)) * 0.6)
		}
 		
#	 		text(xPo,(max(y)*.92),cex=fontSize,pos=4,substitute(r^2 == k, list(k = round(summary(r)$adj.r.squared,2))))
#			text(xPo,(max(y)*.84),cex=fontSize,pos=4,substitute(x == k, list(k = round(j,2))))
#			text(xPo,(max(y)*.76),cex=fontSize,pos=4,substitute(n == k,list (k = length(x))))
		text(xPo,(min(y) + (max(y)-min(y))*.21),cex=fontSize,pos=4,substitute(r^2 == k, list(k = round(summary(r)$adj.r.squared,2))))
		if (j < 0.01) {
			text(xPo,(min(y) + (max(y)-min(y))*.11),cex=fontSize,pos=4,substitute(x < k, list(k = 0.01)))
		} else {
			text(xPo,(min(y) + (max(y)-min(y))*.11),cex=fontSize,pos=4,substitute(x == k, list(k = round(j,2))))
		}
		text(xPo,(min(y) + (max(y)-min(y))*.02),cex=fontSize,pos=4,substitute(n == k,list (k = length(x))))
	}
	
	# PLOT ERROR BARS
	if (peb) {
		for (i in 1:length(x[,'fit'])) {
			lines(c(x[i,'fit'],x[i,'fit']),c(y[i,'upr'],y[i,'lwr']), lwd=0.5)
			lines(c(x[i,'upr'],x[i,'lwr']),c(y[i,'fit'],y[i,'fit']), lwd=0.5)
		}
		points(x[,'fit'], y[,'fit'], col='black', pch=22, cex=1)
	}
	# PLOT UNITY LINE
	if (pul) {
		lines (c(mMin,mMax),c(mMin,mMax), lwd=2, lty=2)
	}
	
	# PLOT UNCERTAINTY
	if (pun) {
		r.pred <- predict(r, int='p')
		r.var <- median(r.pred[,'upr'] - r.pred[,'lwr'])
#		text((min(x^j)),(max(y)*.8),cex=1.2,pos=4,paste('unc =',round(r.var,0)))

		aa.x <- data.frame(x=sort(x))
		aa.pred <- predict(r, int='p', newdata=aa.x)
		aa.conf <- predict(r, int='c', newdata=aa.x)
		matlines(aa.x^j,aa.pred, lty=c(5,1,1), col=c('black','grey66', 'grey66'))
		matlines(aa.x^j,aa.conf, lty=c(1,1), col=c('grey33', 'grey33'))

	}

	# PLOT SAMPLE NUMBERS
 	if (!missing(pl)) {
#		text(x^j,y,cex=.5,substr(as.character(pl),3,20))
		text(x^j,y,cex=.5,as.character(pl))
	}

	# RETURN PREDICTION INTERVAL	
	if (pun)
		return (r.pred)
	else
		return (j)
}


regressOpim <- function(x, y) {

	save(x, y, file ='tmp.Rimage')

	# OPTIMIZE EXPONENT FIT TO DATA
	optPara <- optimize(optimExpX, interval = c(0,6))
	j <- optPara$minimum

	return(j)
}


regressBackstrip  <- function(dates, taxon, colAA, namAA) {


}
