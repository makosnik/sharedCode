library(outliers)

###################################################
#	FORMULAE PROVIDED BY QUAN HUA - 2011.OCT.01
#
#	Bevington PR, Robinson DK. 1992. 
#	Data Reduction and Error Analysis for the Physical Sciences. 
#	Second edition. New York: McGraw-Hill Inc.
###################################################


##	CALCULATE WEIGHTED MEAN
##	-	values to be averaged (x)
##	-	weights (w) should be 1/sigma^2
##	-	length of w and x must be the same.
###################################################
mak.wt.mean <- function (x,w) {
	return(sum(w*x)/sum(w))	
}


##	CALCULATE WEIGHTED STANDARD DEVIATION OF MEAN
##	-	values to be averaged (x)
##	-	weights (w) should be 1/sigma^2
##	-	length of w and x must be the same.
###################################################
mak.wt.sdev <- function (x,w) {
	
	xbar <- mak.wt.mean(x,w)
	xdif <- (x-xbar)^2

	return(sqrt((sum(xdif * w) / sum(w)) * (length(x)/(length(x)-1))))
}


##	CALCULATE ERROR OF THE WEIGHTED MEAN
##	-	weights (w) should be 1/sigma^2
###################################################
mak.wt.errorOfMean <- function (w) {
	return(sqrt(1/sum(w)))	
}


##	CALCULATE DISTANCE BETWEEN UPPER AND LOWER CONFIDENCE DATES
###################################################
mak.intervalSpan <- function (ageLwr,ageUpr,aminos) {
	
	ageUncs <- vector(mode='logical', length(aminos))
	for (amino in 1:length(aminos)) {
		ageUncs[amino] <- paste(aminos[amino],'AgeUnc',sep='_')
	}
	
	##	DETERMINE IF THE LOWEST VALUE IS NEGATIVE TO PROPERLY CALACULATE RANGE
	big <- 0
	lowestValue <- (min(ageLwr, na.rm=TRUE))
	if (lowestValue < 0)
		big <- abs(lowestValue)
	
	##	DEFINE & CALCULATE ALL THE AGE RANGES
	ageUncert <- matrix(nrow=nrow(ageLwr),ncol=length(aminos))
	colnames(ageUncert) <- ageUncs

	for (amino in 1:length(aminos)) {
		ageUncert[,ageUncs[amino]] <- unlist(ageUpr[amino]+big) - unlist(ageLwr[amino]+big)
	}

	return(ageUncert)
}


##	WRAPPER TO CALCULATE ALL OF THE VARIOUS WEIGHTED STATS
##	-	DOES NOT ADD THE NEEDED COLUMNS
##	-	RETURNS THE SAME COLUMNS AS SENT
###################################################
mak.wt.calcs <- function(dataSet,aminos) {
	
	##	DEFINE ALL THE AGE COLUMN NAMES
	ageFits <- vector(mode='logical', length(aminos))
	ageLwrs <- vector(mode='logical', length(aminos))
	ageUprs <- vector(mode='logical', length(aminos))

	for (amino in 1:length(aminos)) {
		ageFits[amino] <- paste(aminos[amino],'AgeFit',sep='_')
		ageLwrs[amino] <- paste(aminos[amino],'AgeLwr',sep='_')
		ageUprs[amino] <- paste(aminos[amino],'AgeUpr',sep='_')
	}

	## GET AGE CONFIDENCE INTERVAL WIDTH
	ageUncert <- mak.intervalSpan(dataSet[,ageLwrs], dataSet[,ageUprs], aminos)
	## CANNOT ALLOW AN UNCERTAINTY OF 0, SET MINIMUM UNCERTAINTY TO 1 YEAR
	ageUncert[(ageUncert[,]==0)] <- 1

	## CALCULATE WEIGHTED MEAN, SD, EM, AND UNCERTAINTY
	for (i in 1:nrow(dataSet)) {
		dataSet[i,'wmAll'] <- round(mak.wt.mean(unlist(dataSet[i,ageFits]),unlist(1/ageUncert[i,]^2)),0)		
		dataSet[i,'sdAll'] <- round(mak.wt.sdev(unlist(dataSet[i,ageFits]),unlist(1/ageUncert[i,]^2)),0)
		dataSet[i,'emAll'] <- round(mak.wt.errorOfMean(unlist(1/ageUncert[i,]^2)),0)		
		dataSet[i,'wuAll'] <- dataSet[i,'sdAll']

		if (dataSet[i,'emAll'] > dataSet[i,'sdAll']) {
			dataSet[i,'wuAll'] <- dataSet[i,'emAll']
		}
	}
	return(dataSet)
}


##	WRAPPER TO CALCULATE ALL OF THE VARIOUS WEIGHTED STATS
##	-	ADDS THE NEEDED COLUMNS
##	-	ALSO RUNS GRUBB'S TESTS
###################################################
mak.wrapper.wt.calcs <- function(dataSet,aminos) {

	##	DEFINE ALL THE AGE COLUMN NAMES
	ageFits <- vector(mode='logical', length(aminos))
	ageLwrs <- vector(mode='logical', length(aminos))
	ageUprs <- vector(mode='logical', length(aminos))
	ageUncs <- vector(mode='logical', length(aminos))

	for (amino in 1:length(aminos)) {
		ageFits[amino] <- paste(aminos[amino],'AgeFit',sep='_')
		ageLwrs[amino] <- paste(aminos[amino],'AgeLwr',sep='_')
		ageUprs[amino] <- paste(aminos[amino],'AgeUpr',sep='_')
		ageUncs[amino] <- paste(aminos[amino],'AgeUnc',sep='_')
	}

	##	SET UP VECTORS TO STORE WEIGHTED ANSWERS
	wmAll <- vector(mode='logical',nrow(dataSet))
	wuAll <- vector(mode='logical',nrow(dataSet))
	sdAll <- vector(mode='logical',nrow(dataSet))
	emAll <- vector(mode='logical',nrow(dataSet))
	dataSet <- cbind(dataSet,wmAll,wuAll,sdAll,emAll)

	##	CALCULATE WEIGHTED MEANS
	dataSet <- mak.wt.calcs(dataSet,aminos)

	##	SET UP VECTORS TO STORE GRUBB'S ANSWERS
	g10.p <- vector(mode='logical',nrow(dataSet))
	g10.h <- vector(mode='logical',nrow(dataSet))
	g10.wm <- vector(mode='logical',nrow(dataSet))
	g10.wu <- vector(mode='logical',nrow(dataSet))
	g10.aa <- vector(mode='logical',nrow(dataSet))

	pCrit <- (0.05/nrow(dataSet))
#	pCrit <- 0.05
	
	##	RUN GRUBB'S TESTS
	for (i in 1:nrow(dataSet)) {

		## REQUIRE 6 OR MORE TO RUN TEST
		if (length(aminos) > 5) {
						
			##	GRUBB'S TESTS
			g10 <- grubbs.test(unlist(dataSet[i,ageFits]), type=10, two.sided=TRUE)
			g10.p[i] <- g10$p.value			
			g10.h[i] <- g10$alternative			
									
			if (g10$p.value < pCrit) {
				g10.aa[i] <- substr(unlist(dimnames(as.array(g10$statistic[1]))),3,5)
				dd <- mak.wt.calcs(dataSet[i,],setdiff(aminos,c(g10.aa[i])))
				g10.wm[i] <- dd[1,'wmAll']
				g10.wu[i] <- dd[1,'wuAll']
			}
		}
	}
	
	dataSet <- cbind(dataSet,g10.p,g10.h,g10.aa,g10.wm,g10.wu)
	
	return(dataSet)
}

###################################################
mak.wrapper.plot.ages <- function(dataSet, aminos, ymin, ymax) {

	alab <- 0.75
	i <- 1

	if (missing(aminos)) {
		aminos <- c('Asp','Glu','Ala','Val','Phe','Leu','A.I')
	}

	ageFit <- vector(mode='logical', length(aminos))
	ageUpr <- vector(mode='logical', length(aminos))
	ageLwr <- vector(mode='logical', length(aminos))
	
	for (amino in 1:length(aminos)) {
		ageFit[amino] <- paste(aminos[amino],'AgeFit',sep='_')
		ageLwr[amino] <- paste(aminos[amino],'AgeLwr',sep='_')
		ageUpr[amino] <- paste(aminos[amino],'AgeUpr',sep='_')
	}


	for (i in 1:nrow(dataSet)) {
		
		if (missing(ymax))
			ymax <- ceiling((max(dataSet[i,c(ageUpr,'X2sOld')], na.rm=TRUE)*1.03)/100)*100
		if (missing(ymin))
			ymin <- floor((min(dataSet[i,c(ageLwr,'X2sYng')], na.rm=TRUE)*0.97)/100)*100
	
		##	PLOT AXES
		plot(1:1,type='n',ylim=c(ymin,ymax),xlim=c(0,13), xaxt='n', ylab='Age (years)', xlab='',xaxs='i',yaxs='i')
		title(paste(dataSet[i,'material'],', specimen #',dataSet[i,'specimen_no'],sep=''),font.main=1)
		
		
		##	PLOT AMINO ACIDS
		a <- 0
		for (amino in 1:length(aminos)) {
			a <- a+1
			points(a,dataSet[i, ageFit[amino]], pch=19)
			lines(c(a,a), c(dataSet[i, ageLwr[amino]],dataSet[i, ageUpr[amino]]))
			text(a+0.3,dataSet[i, ageLwr[amino]], aminos[amino], cex= alab)
		}
		
		##	PLOT CARBON-14 AGE
		a <- a+1
		points(a,dataSet[i,'X2sMedian'], pch=19)
		lines(c(a,a), c(dataSet[i,'X2sYng'],dataSet[i,'X2sOld']))
		text(a+0.3, dataSet[i,'X2sYng'], '14C', cex=alab, pos=1)
		
		##	DIVIDER LINE
		a <- a+1
		lines(c(a,a),c(ymin,ymax))
	
		##	PLOT WEIGHTED MEAN
		a <- a+1
		points(a,dataSet[i,'wmAll'], pch=19)
		lines(c(a,a), c(dataSet[i,'wmAll']+ dataSet[i,'wuAll'], dataSet[i,'wmAll']-dataSet[i,'wuAll']))
		text(a+0.3,(dataSet[i,'wmAll']-dataSet[i,'wuAll']), 'wmean', cex=alab, pos=1)
	
		##	PLOT ASP V. GLU ONLY
		a <- a+1
		dd <- mak.wt.calcs(dataSet[i,],c('Asp','Glu'))
		points(a,dd[i,'wmAll'], pch=19)
		lines(c(a,a), c(dd[i,'wmAll']+dd[i,'wuAll'],dd[i,'wmAll']-dd[i,'wuAll']))
		text(a+0.3,(dd[i,'wmAll']-dd[i,'wuAll']), 'Asp&Glu', cex=alab, pos=1)
	
		##	PLOT -GLU FOR LILOA
		if (dataSet[i,'material']=='Liloa') {
			a <- a+1
			dd <- mak.wt.calcs(dataSet[i,],c('Asp','Ala','Val','Phe','Leu','A.I'))
			points(a,dd[i,'wmAll'], pch=19)
			lines(c(a,a), c(dd[i,'wmAll']+dd[i,'wuAll'],dd[i,'wmAll']-dd[i,'wuAll']))
			text(a+0.3,(dd[i,'wmAll']-dd[i,'wuAll']), '-Glu', cex=alab, pos=1)
		}
	}
}


###################################################
mak.wrapper.varationPlot <- function (dataSet,fileName,xCol,yCol,TAXA,maxSD,maxCV) {
	
	print(paste("writing eps file: ", fileName, sep=''))
	
	postscript(fileName)
	par(mar=c(4,4,2,2), las=0, mgp=c(2,1,0))
	plotSize <- 7
	rows <- 3
	cols <- 2
	nplots <- rows*cols
	layout(matrix(seq(nplots), rows, cols, byrow=TRUE), width=rep(lcm(plotSize),cols), height=rep(lcm(plotSize),rows), respect=TRUE)
	plotLetter <- 1
		
	if (missing(maxSD))
		maxSD <- max(dataSet[,yCol])
	if (missing(maxCV))
		maxCV <- max(abs(dataSet[,yCol]/dataSet[,xCol]),na.rm=TRUE)
	xmax <- max(dataSet[,xCol])

	TAXA <- levels(as.factor(dataSet[,'material']))
	
	for (taxon in TAXA) {
		dataTaxon <- dataSet[(dataSet[,'material']==taxon),]
		mak.varationPlot(dataTaxon[,xCol], dataTaxon[,yCol], plotLetter, taxon, ymax=maxSD)
	
	}
	
	dev.off()
}


###################################################
mak.varationPlot <- function (xData, yData, plotLetter, taxon, ymax) {
	
	xmax <- max(xData)
	if (missing(ymax))
		ymax <- max(yData)

	plot(xData, yData, main='', ylab='uncertainty of the mean (years)',xlab='mean age (years)', pch=21, cex=1, col='white', bg='black', xlim=c(0,xmax),ylim=c(0,ymax))

	title (main=paste(LETTERS[plotLetter],'.',sep=''),adj=0)
	title (main=taxon, adj=0.5, font.main=3)

	# Y lines
	for (y in c(0.05,0.1,0.2,0.4)) {
		x <- seq((1000),(xmax*1.3),length.out=2)
		lines(x,x*y, lty=3)
		z <- c(-100,1000)
		lines(z,rep((1000*y),2),lty=3)
	}
	lines(c(xmax*1.3,-100),c(0,0),lty=3)

	for (y in c(0.05,0.1,0.2,0.4)) {
		x <- seq((1000),(xmax*1.3),length.out=2)
		lines(x,x*y, lty=3)
		z <- c(-100,1000)
		lines(z,rep((1000*y),2),lty=3)
	}
	
}


###################################################
mak.regressPlot <- function (xData, yData, plotLetter, taxon) {
	
	regressPlot(xData,as.matrix(yData),xl='Glu age (years)',yl='Asp age (years)',pul=TRUE, peb=TRUE, pr2=TRUE)
	for (i in 1:length(xData[,'fit'])) {
		lines(c(xData[i,'fit'], xData[i,'fit']),c(yData[i,'upr'],yData[i,'lwr']), lwd=1.0)
		lines(c(xData[i,'upr'], xData[i,'lwr']),c(yData[i,'fit'],yData[i,'fit']), lwd=1.0)
	}
	title (main=paste(LETTERS[plotLetter],'.',sep=''),adj=0)
	title (main=taxon, adj=0.5, font.main=3)
	
}

