#########################################################################################
##	ADD AAR CONCENTRATIONS AND RATIOS COLUMNS TO DATA TABLE
#########################################################################################
addaarConcColumns <- function(areaRaw) {


	###################################################
	#	CALCULATE - D/L RATIOS
	###################################################

	areaRaw <- cbind(areaRaw,cD_Asp, cL_Asp, cD_Glu, cL_Glu, cD_Ser, cL_Ser, cD_Ala, cL_Ala, cD_Val, cL_Val, cD_Phe, cL_Phe, cD_Leu, cL_Leu, cD_AI, cL_AI)
	areaRaw <- cbind(areaRaw, Asp_DL, Glu_DL, Ser_DL, Ala_DL, Val_DL, Phe_DL, Leu_DL, AI_DL)
	
	return(areaRaw)
}

#########################################################################################
##	ADD AAR CONCENTRATIONS AND RATIOS COLUMNS TO DATA TABLE
#########################################################################################
addaarDLColumns <- function(areaRaw) {

	###################################################
	#	CALCULATE - D & L CONCENTRATIONS
	###################################################
	cD_Asp <- (areaRaw[,'D_Asp']/areaRaw[,'L_hArg'])*200
	cL_Asp <- (areaRaw[,'L_Asp']/areaRaw[,'L_hArg'])*200
	cD_Glu <- (areaRaw[,'D_Glu']/areaRaw[,'L_hArg'])*200
	cL_Glu <- (areaRaw[,'L_Glu']/areaRaw[,'L_hArg'])*200
	cD_Ser <- (areaRaw[,'D_Ser']/areaRaw[,'L_hArg'])*200
	cL_Ser <- (areaRaw[,'L_Ser']/areaRaw[,'L_hArg'])*200
	cD_Ala <- (areaRaw[,'D_Ala']/areaRaw[,'L_hArg'])*200
	cL_Ala <- (areaRaw[,'L_Ala']/areaRaw[,'L_hArg'])*200
	cD_Val <- (areaRaw[,'D_Val']/areaRaw[,'L_hArg'])*200
	cL_Val <- (areaRaw[,'L_Val']/areaRaw[,'L_hArg'])*200
	cD_Phe <- (areaRaw[,'D_Phe']/areaRaw[,'L_hArg'])*200
	cL_Phe <- (areaRaw[,'L_Phe']/areaRaw[,'L_hArg'])*200
	cD_Leu <- (areaRaw[,'D_Leu']/areaRaw[,'L_hArg'])*200
	cL_Leu <- (areaRaw[,'L_Leu']/areaRaw[,'L_hArg'])*200
	cD_Ile  <- (areaRaw[,'D_AIle']/areaRaw[,'L_hArg'])*200
	cL_Ile  <- (areaRaw[,'L_Ile']/areaRaw[,'L_hArg'])*200

	areaRaw <- cbind(areaRaw,cD_Asp, cL_Asp, cD_Glu, cL_Glu, cD_Ser, cL_Ser, cD_Ala, cL_Ala, cD_Val, cL_Val, cD_Phe, cL_Phe, cD_Leu, cL_Leu, cD_Ile, cL_Ile)

	#	CALCULATE - TOTAL CONCENTRATIONS
	###################################################
	cT_Asp <- areaRaw[,'cD_Asp']+areaRaw[,'cL_Asp']
	cT_Glu <- areaRaw[,'cD_Glu']+areaRaw[,'cL_Glu']
	cT_Ser <- areaRaw[,'cD_Ser']+areaRaw[,'cL_Ser']
	cT_Ala <- areaRaw[,'cD_Ala']+areaRaw[,'cL_Ala']
	cT_Val <- areaRaw[,'cD_Val']+areaRaw[,'cL_Val']
	cT_Phe <- areaRaw[,'cD_Phe']+areaRaw[,'cL_Phe']
	cT_Leu <- areaRaw[,'cD_Leu']+areaRaw[,'cL_Leu']
	cT_Ile <- areaRaw[,'cD_Ile']+areaRaw[,'cL_Ile']

areaRaw <- cbind(areaRaw, cT_Asp, cT_Glu, cT_Ser, cT_Ala, cT_Val, cT_Phe, cT_Leu, cT_Ile)

	###################################################
	#	CALCULATE - D/L RATIOS
	###################################################
	Asp_DL <- areaRaw[,'D_Asp']/areaRaw[,'L_Asp']
	Glu_DL <- areaRaw[,'D_Glu']/areaRaw[,'L_Glu']
	Ser_DL <- areaRaw[,'D_Ser']/areaRaw[,'L_Ser']
	Ala_DL <- areaRaw[,'D_Ala']/areaRaw[,'L_Ala']
	Val_DL <- areaRaw[,'D_Val']/areaRaw[,'L_Val']
	Phe_DL <- areaRaw[,'D_Phe']/areaRaw[,'L_Phe']
	Leu_DL <- areaRaw[,'D_Leu']/areaRaw[,'L_Leu']
	Ile_DL <- areaRaw[,'D_AIle']/areaRaw[,'L_Ile']

	areaRaw <- cbind(areaRaw, Asp_DL, Glu_DL, Ser_DL, Ala_DL, Val_DL, Phe_DL, Leu_DL, Ile_DL)
	
	return(areaRaw)
}


#########################################################################################
##	AMINO ACID D/L CORRELATIONS
#########################################################################################
plotAARCrossCorrelations <- function(dataToUse, iniA3, fileName, highlightPoints) {
	
	if (!is.na(fileName))
		pdf(paste(PATH.OUT,fileName,sep='/'), width=pageWidthTwo, height=pageWidthTwo)
	par(mfcol=c(length(iniA3),length(iniA3)), cex=0.5, mar=c(1,1,1,1), mgp=c(1,0.5,0))
	
	resid <- matrix(nrow=nrow(dataToUse),ncol=length(iniA3)*length(iniA3))
	colname <- list()
	i <- 0
	for (x in iniA3) {
		for (y in iniA3) {
			i<-i+1
			colname[i] <- paste(x,y,sep='')
		}
	}
	colnames(resid) <- colname

	for (x in iniA3) {
		for (y in iniA3) {
		
			## GET X AND Y DATA
			xData <- log(dataToUse[,paste(x,'DL',sep='_')])
			yData <- log(dataToUse[,paste(y,'DL',sep='_')])
			xData <- (dataToUse[,paste(x,'DL',sep='_')])
			yData <- (dataToUse[,paste(y,'DL',sep='_')])
		
			## IF X = Y THEN PLOT A LABEL BOX RATHER THAN X v. Y
			if (x==y) {
			
				## EMPTY PLOT BOX
				plot(xData, yData, xaxt='n', yaxt='n', xlab='', ylab='', lwd=0.2, type='n')
			
				# SKIP SIDES 2 & 3 FOR UPPER LEFT CORNER
				if (x != iniA3[1]) {
					axis(2,at=NULL)
					axis(3,at=NULL)
				}
				# SKIP SIDES 1 & 4 FOR THE LOWER RIGHT CORNER
				if (x != iniA3[length(iniA3)]) {
					axis(1,at=NULL)
					axis(4,at=NULL)
				}
			
				# LABEL THE AMINMO ACID IN THE CENTER OF THE PLOT SPACE
				cenPt <- ((max(xData, na.rm=TRUE)-min(xData, na.rm=TRUE))/2)+min(xData, na.rm=TRUE)
				if (x !='AI')
					text(cenPt, cenPt, paste(x,'D/L'), cex=2, pos=3)
				else
					text(cenPt, cenPt, 'A/I', cex=2, pos=3)			
				text(cenPt, cenPt, taxon, cex=1.5, font=3, pos=1)

			## PLOT X v. Y
			} else {
			
				plot(xData, yData, xaxt='n', yaxt='n', xlab='', ylab='', pch=21, col='white', bg='grey', lwd=0.2)
				text(xData, yData, dataToUse$specimen_no, cex=0.25, pos=4)
				
	#			if (!is.na(highlightPoints))
				for (h in highlightPoints) {
					points(xData[(dataToUse$specimen_no == h)], yData[(dataToUse$specimen_no == h)], pch=21, col='blue',cex=1.2)
					text(xData[(dataToUse$specimen_no == h)], yData[(dataToUse$specimen_no == h)],dataToUse[(dataToUse$specimen_no == h),'specimen_no'], col='blue', pos=2)
				}
								
				lines(c(0,1),c(0,1), lty=2, lwd=0.5)
			
				# CALCULATE CORRELATIONS
				a <- cor.test(xData, yData, method='spearman')
				b <- lm(yData ~ xData)
				resid[, paste(x,y,sep='')] <- b$residuals
			
				# PLOT CORRELATIONS
				abline(b)
				xPad <- ((max(xData, na.rm=TRUE)-min(xData, na.rm=TRUE))/20)
				yPad <- ((max(yData, na.rm=TRUE)-min(yData, na.rm=TRUE))/30)
				text(max(xData, na.rm=TRUE)+xPad,min(yData, na.rm=TRUE)+yPad, substitute(rho==j,list(j=round(a$estimate,2))), pos=2)
				text(min(xData, na.rm=TRUE)-xPad,max(yData, na.rm=TRUE)-yPad, substitute(r^2==k, list(k=round(summary(b)$r.squared,2))), pos=4)
					
			}
		}
	}
	if (!is.na(fileName))
		dev.off()
	return(resid)
}

#########################################################################################
#	PLOT AMINO ACID D/L RESIDUALS
#		- need to fix this to:
#			- take arbretary list of Amino Acids
#			- dynamic y-axis
#########################################################################################
plotAARResiduals <- function(dataOne, resid, iniA4, fileName) {

	quants <- matrix(nrow=nrow(dataOne),ncol=5)
	for (i in 1:nrow(resid)) {
		quants[i,] <- quantile((resid[i,]), c(0.025,0.25,0.5,0.75,0.975), na.rm=TRUE)
	}

	if (!is.na(fileName))
		pdf(paste(PATH.OUT,fileName,sep='/'), width=pageWidthTwo, height=pageWidthTwo)
	
	par(mfcol=c(1,1), mar=c(8,2,1,1))
	
	plot(c(1,nrow(resid)),c(-0.165,max(resid, na.rm=TRUE)), type='n', xlab='',xaxt='n')
	axis(1,at=seq(1,nrow(resid), by=1),labels=paste(dataOne[,'ual_num'],dataOne[,'specimen_no'],sep=' - '), las=3)
	lines(c(-1,nrow(resid))+2,c(0,0))
	segments(x0=seq(1,nrow(resid), by=1),y0=quants[,1],y1=quants[,5], col='grey')
	segments(x0=seq(1,nrow(resid), by=1),y0=quants[,2],y1=quants[,4], lwd=5, col='grey')
	for (x in iniA4) {
		for (y in iniA4) {
			points(seq(1,nrow(resid), by=1),(resid[,paste(x,y,sep='')]), pch=19, cex=0.5)
			text(seq(1,nrow(resid), by=1),(resid[,paste(x,y,sep='')]), paste(x,y,sep=''), cex=0.5, pos=4)
	}}
	
	text(-1,-0.14,'Asp', cex=0.6)
	text(seq(1,nrow(resid), by=1),-0.14, round(dataOne[,'Asp_DL'],2), cex=0.4)
	text(-1,-0.145,'Glu', cex=0.6)
	text(seq(1,nrow(resid), by=1),-0.145, round(dataOne[,'Glu_DL'],2), cex=0.4)
	text(-1,-0.15,'Ala', cex=0.6)
	text(seq(1,nrow(resid), by=1),-0.15, round(dataOne[,'Ala_DL'],2), cex=0.4)
	text(-1,-0.155,'Val', cex=0.6)
	text(seq(1,nrow(resid), by=1),-0.155, round(dataOne[,'Val_DL'],2), cex=0.4)
	text(-1,-0.16,'Phe', cex=0.6)
	text(seq(1,nrow(resid), by=1),-0.16, round(dataOne[,'Phe_DL'],2), cex=0.4)
	
	if (!is.na(fileName))
		dev.off()
}


#########################################################################################
#	PLOT AMINO ACID D/L REPLICATES
#########################################################################################
plotAARReplicates <- function(dataOne, iniA3, fileName) {

	specimens <- unique(dataOne[,'specimen_no'])
	dataReps <- dataOne[1,]
	
	if (length(dataOne[,'specimen_no']) == length(unique(dataOne[,'specimen_no'])))
		return(-1)
	
	for (s in specimens) {
		tmp <- dataOne[(dataOne[,'specimen_no']==s),]
		if (nrow(tmp) > 1)
			dataReps <- rbind(dataReps,tmp)
	}
	dataReps <- dataReps[-1,]
	
	if (!is.na(fileName))
		pdf(paste(PATH.OUT,fileName,sep='/'), width=pageWidthTwo, height=pageWidthTwo)
	
	par(mfcol=c(length(iniA3),length(iniA3)), cex=0.5, mar=c(1,1,1,1), mgp=c(1,0.5,0))

	diagonalLabel <- '(replicate specimens)'


	for (x in iniA3) {
		for (y in iniA3) {
	
			## GET X AND Y DATA
			xData <- dataToUse[,paste(x,'DL',sep='_')]
			yData <- dataToUse[,paste(y,'DL',sep='_')]
		
			## IF X = Y THEN PLOT A LABEL BOX RATHER THAN X v. Y
			if (x==y) {

				## EMPTY PLOT BOX
				plot(xData, yData, xaxt='n', yaxt='n', xlab='', ylab='', lwd=0.2, type='n')
			
				# SKIP SIDES 2 & 3 FOR UPPER LEFT CORNER
				if (x != iniA3[1]) {
					axis(2,at=NULL)
					axis(3,at=NULL)
				}
				# SKIP SIDES 1 & 4 FOR THE LOWER RIGHT CORNER
				if (x != iniA3[length(iniA3)]) {
					axis(1,at=NULL)
					axis(4,at=NULL)
				}
			
				# LABEL THE AMINMO ACID IN THE CENTER OF THE PLOT SPACE
				cenPt <- ((max(xData, na.rm=TRUE)-min(xData, na.rm=TRUE))/2)+min(xData, na.rm=TRUE)
				if (x !='AI')
					text(cenPt, cenPt, paste(x,'D/L'), cex=2, pos=3)
				else
					text(cenPt, cenPt, 'A/I', cex=2, pos=3)			
				text(cenPt, cenPt, taxon, cex=1.5, font=3, pos=1)

			## PLOT X v. Y
			} else {
				
				plot(xData, yData, xaxt='n', yaxt='n', xlab='', ylab='', pch=21, col='white', bg='grey', lwd=0.2)
				lines(c(0,1),c(0,1), lty=2, lwd=0.5)
			
				# CALCULATE CORRELATIONS
				a <- cor.test(xData, yData, method='spearman')
				b <- lm(yData ~ xData)

				# PLOT CORRELATIONS
				abline(b)
				xPad <- ((max(xData, na.rm=TRUE)-min(xData, na.rm=TRUE))/20)
				yPad <- ((max(yData, na.rm=TRUE)-min(yData, na.rm=TRUE))/20)
				text(max(xData, na.rm=TRUE)+xPad,min(yData, na.rm=TRUE)+yPad, substitute(rho==j,list(j=round(a$estimate,2))), pos=2)
				text(min(xData, na.rm=TRUE)-xPad,max(yData, na.rm=TRUE)-yPad, substitute(r^2==k, list(k=round(summary(b)$r.squared,2))), pos=4)
							
				specs <- unique(dataReps[,'specimen_no'])
				for (s in specs) {
					dataSpec <- dataReps[(dataReps[,'specimen_no']==s),]
#					points(dataReps[,paste(x,'DL',sep='_')], dataReps[,paste(y,'DL',sep='_')], pch=21, col='black')
					text(dataSpec[1,paste(x,'DL',sep='_')], dataSpec[1,paste(y,'DL',sep='_')], dataSpec[1,'specimen_no'], pos=4, cex=0.5)
					if (nrow(dataSpec) > 1) {
						for (r in 2:(nrow(dataSpec))) {
							segments(dataSpec[r,paste(x,'DL',sep='_')], dataSpec[r,paste(y,'DL',sep='_')],dataSpec[(r-1),paste(x,'DL',sep='_')], dataSpec[(r-1),paste(y,'DL',sep='_')], lwd=0.5)
						}
					}
					if (nrow(dataSpec) > 2) {
						segments(dataSpec[1,paste(x,'DL',sep='_')], dataSpec[1,paste(y,'DL',sep='_')],dataSpec[(nrow(dataSpec)),paste(x,'DL',sep='_')], dataSpec[(nrow(dataSpec)),paste(y,'DL',sep='_')], lwd=0.5)
					}
				}
			}
		}
	}

	
	if (!is.na(fileName))
		dev.off()
}


#########################################################################################
##	PLOT AAR FRACTIONS
#########################################################################################
plotAARfractions <- function (wholeData, aminoList, xFraction, yFraction, fileName) {

	xData <- wholeData[(wholeData[,'analysis_fraction'] == xFraction),]
	yData <- wholeData[(wholeData[,'analysis_fraction'] == yFraction),]
	
	yLab <- yFraction
	xLab <- xFraction
	
	if (xLab == '')
		xLab <- 'THAA'
	if (yLab == '')
		yLab <- 'THAA'

	if (yFraction == 'B-FAA')
		yLab <- 'Intra-crystalline FAA'
	if (yFraction == 'B-THAA')
		yLab <- 'Intra-crystalline THAA'
	if (xFraction == 'B-FAA')
		xLab <- 'Intra-crystalline FAA'
	if (xFraction == 'B-THAA')
		xLab <- 'Intra-crystalline THAA'

	xLab <- paste(xFraction,'D/L')
	yLab <- paste(yFraction,'D/L')
		
	jointSpecs <- intersect(xData$specimen_no, yData$specimen_no)
	
	if (!is.vector(jointSpecs))
		return(-1)
	
	xData <- xData[(xData[,'specimen_no'] %in% jointSpecs),]
	yData <- yData[(yData[,'specimen_no'] %in% jointSpecs),]

	xMin <- min(xData[,paste(aminoList,'DL',sep='_')], na.rm=TRUE)
	xMax <- max(xData[,paste(aminoList,'DL',sep='_')], na.rm=TRUE)
	yMin <- min(yData[,paste(aminoList,'DL',sep='_')], na.rm=TRUE)
	yMax <- max(yData[,paste(aminoList,'DL',sep='_')], na.rm=TRUE)
	dMin <- floor(min(xMin,yMin)*10)/10
	dMax <- ceiling(max(xMax,yMax)*10)/10

	if (!is.na(fileName))
		pdf(paste(PATH.OUT,fileName,sep='/'), width=pageWidthOne, height=pageHeight)

#	par(mfcol=c(ceiling(length(aminoList)/2),2), cex=1, mar=c(4,4,1,1), mgp=c(2,0.7,0), cex.lab=1.25, las=1)
	par(mfcol=c(3,1), cex=1, mar=c(4,4,1,1), mgp=c(2,0.7,0), oma=c(1,1,1,1),cex.lab=1.25, las=1)
	par(mfcol=c(3,1), cex=1, mar=c(3.1,3.5,0,0), mgp=c(2.1,0.7,0), oma=c(0.5,0.5,1,1),cex.lab=1.25, las=1)

	yLT <- aminoList[seq(1, ceiling(length(aminoList)/2), by=1)]
	xLT <- aminoList[c(ceiling(length(aminoList)/2), length(aminoList))]
	xLT <- aminoList[length(aminoList)]
	
	for (aa in aminoList) {
	
		DL <- paste(aa,'DL',sep='_')
		
#		ylab <- ''
#		if (aa %in% yLT)
			ylab <- yLab
		xlab <- ''
		if (aa %in% xLT)
			xlab <- xLab
		
#		if (aa == 'Asp')
#			par(mar=c(3,3,0,0))
#		if (aa == 'Ala')
#			par(mar=c(2,2.5,1,0.5))
#		if (aa == 'Glu')
#			par(mar=c(3,3,0,0))
#		if (aa == 'Phe')
#			par(mar=c(3,2.5,0,0.5))
		
		
		plot(1:1, type='n', ylim=c(dMin,dMax),xlim=c(dMin,dMax), ylab=ylab, xlab=xlab, xaxs='i', yaxs='i')
		segments(-1,-1,1.5,1.5, lty=1)
		text(dMin, dMax*0.9 ,paste(aa,'D/L'), pos=4, cex=1.25)

		xyPairs <- c(0,0)

		for (s in jointSpecs) {
			xSpec <- xData[(xData[,'specimen_no']==s),]
			ySpec <- yData[(yData[,'specimen_no']==s),]
#			text(median(xSpec[,DL], na.rm=TRUE), median(ySpec[,DL], na.rm=TRUE), xSpec[1,'specimen_no'], pos=1, cex=0.5)
			
			xSq <- quantile(xSpec[,DL], qSD)
			ySq <- quantile(ySpec[,DL], qSD)
			
			points(xSq[3], ySq[3], pch=19, cex=0.5)
			segments(xSq[3], ySq[1], xSq[3], ySq[5], lwd=0.5)
			segments(xSq[1], ySq[3], xSq[5], ySq[3], lwd=0.5)
			
			xyPairs <- rbind(xyPairs,c(xSq[3], ySq[3]))
		}

		xyPairs <- xyPairs[-1,]
		
		lm.p <- lm(xyPairs[,2] ~ xyPairs[,1])
		abline(lm.p, lty=2)
		
		
		rs <- format(round(summary(lm.p)$r.squared,3),nsmall=3)
		text(dMax,0.22, substitute(r^2==k, list(k=rs)), pos=2, cex=0.8)
		psm <- format(round(summary(lm.p)$coefficients[2,1],2),nsmall=2)
		pse <- format(round(summary(lm.p)$coefficients[2,2],2),nsmall=2)
		text(dMax,0.14,paste('slope =',psm,'\u00b1',pse), pos=2, cex=0.8)
		pim <- format(round(summary(lm.p)$coefficients[1,1],2),nsmall=2)
		pie <- format(round(summary(lm.p)$coefficients[1,2],2),nsmall=2)
		text(dMax,0.06,paste('intercept =',pim,'\u00b1',pie), pos=2, cex=0.8)	
	}
	
	if (!is.na(fileName))
		dev.off()
	return()
}


#########################################################################################
##	PLOT AAR FRACTIONS
#########################################################################################
plotAARfractionLossVdl <- function (wholeData, aminoList, xFraction, yFraction, fileName) {

	xData <- wholeData[(wholeData[,'analysis_fraction'] == xFraction),]
	yData <- wholeData[(wholeData[,'analysis_fraction'] == yFraction),]
	
	for (aa in aminoList)
		yData[, paste(aa,'Loss',sep='_')] <- (yData[,paste('cT',aa,sep='_')]/xData[,paste('cT',aa,sep='_')]) * 100
	
	yLab <- yFraction
	xLab <- xFraction
	
#	if (xLab == '')
		xLab <- 'THAA'
#	if (yLab == '')
		yLab <- 'THAA'

	yLab <- paste('( [',yFraction,'] / [', xFraction,'] ) * 100',sep='')
	xLab <- paste(xFraction,'D/L')

		
	jointSpecs <- intersect(xData$specimen_no, yData$specimen_no)
	
	if (!is.vector(jointSpecs))
		return(-1)
	
	xData <- xData[(xData[,'specimen_no'] %in% jointSpecs),]
	yData <- yData[(yData[,'specimen_no'] %in% jointSpecs),]

	xMin <- min(xData[,paste(aminoList,'DL',sep='_')], na.rm=TRUE)
	xMax <- max(xData[,paste(aminoList,'DL',sep='_')], na.rm=TRUE)
	yMin <- min(yData[,paste(aminoList,'Loss',sep='_')], na.rm=TRUE)
	yMax <- max(yData[,paste(aminoList,'Loss',sep='_')], na.rm=TRUE)

	if (!is.na(fileName))
		pdf(paste(PATH.OUT,fileName,sep='/'), width=pageWidthTwo, height=pageWidthTwo)

	par(mfcol=c(ceiling(length(aminoList)/2),2), cex=1, mar=c(4,4,1,1), mgp=c(2,0.7,0), cex.lab=1.25, las=1)
	yLT <- aminoList[seq(1, ceiling(length(aminoList)/2), by=1)]
	xLT <- aminoList[c(ceiling(length(aminoList)/2), length(aminoList))]
	
	for (aa in aminoList) {
			
		ylab <- ''
		if (aa %in% yLT)
			ylab <- yLab
		xlab <- ''
		if (aa %in% xLT)
			xlab <- xLab
		
		if (aa == 'Asp')
			par(mar=c(2,3,1,0))
		if (aa == 'Ala')
			par(mar=c(2,2.5,1,0.5))
		if (aa == 'Glu')
			par(mar=c(3,3,0,0))
		if (aa == 'Phe')
			par(mar=c(3,2.5,0,0.5))
		
		
		plot(1:1, type='n', ylim=c(yMin,yMax),xlim=c(xMin,xMax), ylab=ylab, xlab=xlab)
		text(xMin, yMax*0.9 ,paste(aa), pos=4, cex=1.25)

		xyPairs <- c(0,0)

		for (s in jointSpecs) {
			xSpec <- xData[(xData[,'specimen_no']==s),]
			ySpec <- yData[(yData[,'specimen_no']==s),]
#			text(median(xSpec[,DL], na.rm=TRUE), median(ySpec[,DL], na.rm=TRUE), xSpec[1,'specimen_no'], pos=1, cex=0.5)
			
			xSq <- quantile(xSpec[, paste(aa,'DL',sep='_')], qSD)
			ySq <- quantile(ySpec[, paste(aa,'Loss',sep='_')], qSD)
			
			points(xSq[3], ySq[3], pch=19, cex=0.5)
			segments(xSq[3], ySq[1], xSq[3], ySq[5], lwd=0.5)
			segments(xSq[1], ySq[3], xSq[5], ySq[3], lwd=0.5)
			
			xyPairs <- rbind(xyPairs,c(xSq[3], ySq[3]))
		}

		xyPairs <- xyPairs[-1,]
		
		lm.p <- lm(xyPairs[,2] ~ xyPairs[,1])
		
		pValue <- coefficients(summary(lm.p))[2,4]
		
		text(xMin, yMax*0.8, substitute(p==k, list(k=round(pValue,2))), pos=4)

		if (pValue < 0.05 ) {
			abline(lm.p, lty=2)
			text(xMin, yMax*0.8, substitute(r^2==k, list(k=round(summary(lm.p)$r.squared,2))), pos=4)
			text(xMax, yMax*0.8,paste('slope =',round(summary(lm.p)$coefficients[2,1],0),'\u00b1',round(summary(lm.p)$coefficients[2,2],0)), pos=2)
			text(xMax, yMax*0.7,paste('intercept =',round(summary(lm.p)$coefficients[1,1],0),'\u00b1',round(summary(lm.p)$coefficients[1,2],0)), pos=2)	
		} else {
			quants <- round(quantile(xyPairs[,2], c(0.025,0.5,0.975)))
			medianStr <- paste(quants[2],'% (',quants[1],'-',quants[3],'%)', sep='')
			text(xMin, yMax*0.7, medianStr, pos=4)
		}
	}
	
	if (!is.na(fileName))
		dev.off()
	return()
}


plotAARfractionsSubPlot <- function(aa,yLim,xLim,ylab,xlab,xData,yData,xVal,yVal,jointSpecs) {
	
		plot(1:1, type='n', ylim=c(yLim),xlim=c(xLim), ylab=ylab, xlab=xlab, log='xy')
#		segments(xLim[1],yLim[1],xLim[2],yLim[2], lty=2)
		text(xLim[1], yLim[2]*0.9 ,aa, pos=4, cex=1.25)

		xyPairs <- c(0,0)

		for (s in jointSpecs) {
			xSpec <- xData[(xData[,'specimen_no']==s),]
			ySpec <- yData[(yData[,'specimen_no']==s),]

			points(xSpec[1, xVal], ySpec[1, yVal], pch=19, cex=0.5)
			if (nrow(xSpec) > 1) {
				for (r in 2:(nrow(xSpec))) {
					segments(xSpec[r, xVal], ySpec[1, yVal],xSpec[(r-1), xVal], ySpec[1, yVal], lwd=0.5)
					points(xSpec[r, xVal], ySpec[1, yVal], pch=19, cex=0.5)
				}
			}
			if (nrow(xSpec) > 2) {
				segments(xSpec[1, xVal], ySpec[1, yVal],xSpec[(nrow(xSpec)), xVal], ySpec[(nrow(ySpec)), yVal], lwd=0.5)
			}
			
			xyPairs <- rbind(xyPairs,c(xSpec[1, xVal], ySpec[1, yVal]))
		}
		
		xyPairs <- xyPairs[-1,]
		
		lm.p <- lm(xyPairs[,2] ~ xyPairs[,1])
		abline(lm.p, lty=2)
		
	return()
}


#########################################################################################
##	PLOT AAR FRACTIONS
#########################################################################################
plotAARfractions2 <- function (wholeData, aminoList, xFraction, yFraction1, yFraction2, colPrefix, colSuffix, fileName) {

	xData <- wholeData[(wholeData[,'analysis_fraction'] == xFraction),]
	yData1 <- wholeData[(wholeData[,'analysis_fraction'] == yFraction1),]
	yData2 <- wholeData[(wholeData[,'analysis_fraction'] == yFraction2),]
	
	yLab1 <- yFraction1
	yLab2 <- yFraction2
	xLab <- xFraction
	
	if (xLab == '')
		xLab <- 'THAA'
	if (yLab1 == '')
		yLab1 <- 'THAA'
	if (yLab2 == '')
		yLab2 <- 'THAA'
		
	jointSpecs <- intersect(xData$specimen_no, yData1$specimen_no)
	jointSpecs <- intersect(jointSpecs, yData2$specimen_no)
	
	if (!is.vector(jointSpecs))
		return(-1)
	
	xData <- xData[(xData[,'specimen_no'] %in% jointSpecs),]
	yData1 <- yData1[(yData1[,'specimen_no'] %in% jointSpecs),]
	yData2 <- yData2[(yData2[,'specimen_no'] %in% jointSpecs),]

	xMin <- min(xData[,paste(colPrefix,aminoList,colSuffix,sep='')], na.rm=TRUE)
	xMax <- max(xData[,paste(colPrefix,aminoList,colSuffix,sep='')], na.rm=TRUE)
	yMin1 <- min(yData1[,paste(colPrefix,aminoList,colSuffix,sep='')], na.rm=TRUE)
	yMax1 <- max(yData1[,paste(colPrefix,aminoList,colSuffix,sep='')], na.rm=TRUE)
	yMin2 <- min(yData2[,paste(colPrefix,aminoList,colSuffix,sep='')], na.rm=TRUE)
	yMax2 <- max(yData2[,paste(colPrefix,aminoList,colSuffix,sep='')], na.rm=TRUE)
	dMin <- floor(min(xMin,yMin1,yMin2)*10)/10
	dMax <- ceiling(max(xMax,yMax1,yMax2)*10)/10
	
	if (!is.na(fileName))
		pdf(paste(PATH.OUT,fileName,sep='/'), width=pageHeight, height=pageHeight/2)

	par(mfcol=c(2,(length(aminoList))), cex=1, mar=c(3,3,0.5,0.5), mgp=c(2,0.7,0), cex.lab=1.25)
	yLT <- aminoList[1]
	xLT <- aminoList
	
	for (aa in aminoList) {
	
		cVal <- paste(colPrefix,aa,colSuffix,sep='')

		xMin <- min(xData[, cVal], na.rm=TRUE)
		xMax <- max(xData[, cVal], na.rm=TRUE)
		yMin1 <- min(yData1[, cVal], na.rm=TRUE)
		yMax1 <- max(yData1[, cVal], na.rm=TRUE)
		yMin2 <- min(yData2[, cVal], na.rm=TRUE)
		yMax2 <- max(yData2[, cVal], na.rm=TRUE)
		dMin <- floor(min(xMin,yMin1,yMin2)*10)/10
		dMax <- ceiling(max(xMax,yMax1,yMax2)*10)/10
		dLim <- c(dMin,dMax)
		
		ylab <- ''
		if (aa %in% yLT)
			ylab <- yLab1
		xlab <- ''
		if (aa %in% xLT)
			xlab <- xLab
		
		plotAARfractionsSubPlot(aa, dLim, dLim,ylab,xlab,xData,yData1,cVal,cVal,jointSpecs)
				
		ylab <- ''
		if (aa %in% yLT)
			ylab <- yLab2

		plotAARfractionsSubPlot(aa, dLim, dLim,ylab,xlab,xData,yData2,cVal,cVal,jointSpecs)
		
	}
	if (!is.na(fileName))
		dev.off()
	return()
}


#########################################################################################
##	PLOT AAR FRACTIONS
#########################################################################################
plotAARfractionsVage <- function (wholeData, aminoList, ageCol, yFractions, yPrefix, ySuffix, fileName) {

	xData <- wholeData[!is.na(wholeData[, ageCol]),]
	
#	source(paste('..','makFunctions.R',sep='/'))
	DATA.RAW <- usePreferedCarbonAges(DATA.RAW) 
	
#	wholeData[, paste(yPrefix, aminoList,ySuffix,sep='')] <- log2(wholeData[, paste(yPrefix, aminoList,ySuffix,sep='')])
	xData[,ageCol] <-(xData[,ageCol]) + (2013 - 1950)
	
	
	yDataBF <- wholeData[(wholeData[,'analysis_fraction'] == 'B-FAA'),]
	yDataBF  <- useMedianHPLCrunValue(yDataBF, paste(yPrefix, aminoList,ySuffix,sep=''), 'specimen_no')

	yDataBT <- wholeData[(wholeData[,'analysis_fraction'] == 'B-THAA'),]
	yDataBT  <- useMedianHPLCrunValue(yDataBT, paste(yPrefix, aminoList,ySuffix,sep=''), 'specimen_no')
	
	yDataAT <- wholeData[(wholeData[,'analysis_fraction'] == 'THAA'),]
	yDataAT  <- useMedianHPLCrunValue(yDataAT, paste(yPrefix, aminoList,ySuffix,sep=''), 'specimen_no')
	
	xLab <- 'Age (yr)'
	xLab <- ''
	
	xMin <- min(xData[, ageCol], na.rm=TRUE)
	xMax <- max(xData[, ageCol], na.rm=TRUE)
	xMin <- floor((xMin)/100)*100
	xMax <- ceiling((xMax)/100)*100
	xMin <- 5
	xLim <- c(xMin,xMax)

	if (!is.na(fileName))
		pdf(paste(PATH.OUT,fileName,sep='/'), width=pageHeight, height=pageHeight/3*2)

	par(mfcol=c(length(yFractions),length(aminoList)), cex=1, mar=c(3,3,0.5,0.5), mgp=c(2,0.7,0), cex.lab=1.25)
	
	for (aa in aminoList) {
	
		cVal <- paste(yPrefix,aa,ySuffix,sep='')

		yMin1 <- min(yDataBF[, cVal], na.rm=TRUE)
		yMax1 <- max(yDataBF[, cVal], na.rm=TRUE)
		yMin2 <- min(yDataBT[, cVal], na.rm=TRUE)
		yMax2 <- max(yDataBT[, cVal], na.rm=TRUE)
		yMin3 <- min(yDataAT[, cVal], na.rm=TRUE)
		yMax3 <- max(yDataAT[, cVal], na.rm=TRUE)
		yMin <- floor(min(yMin1, yMin2, yMin3)*100)/100
		yMax <- ceiling(max(yMax1, yMax2, yMin3)*10)/10
		yLim <- c(yMin,yMax)
		
		if (aa != 'Asp') {
			yLab <- ''
		}
		if (aa == 'Asp')
			yLab <- yFractions[1]
		jointSpecs <- intersect(xData$specimen_no, yDataBF$specimen_no)
		plotAARfractionsSubPlot(aa, yLim, xLim, yLab, '', xData, yDataBF, ageCol, cVal, jointSpecs)

		if (aa == 'Asp')
			yLab <- yFractions[2]
		jointSpecs <- intersect(xData$specimen_no, yDataBT$specimen_no)
		plotAARfractionsSubPlot(aa, yLim, xLim, yLab, '', xData, yDataBT, ageCol, cVal, jointSpecs)

		if (aa == 'Asp')
			yLab <- yFractions[3]
		jointSpecs <- intersect(xData$specimen_no, yDataAT$specimen_no)
		plotAARfractionsSubPlot(aa, yLim, xLim, yLab, xLab, xData, yDataAT, ageCol, cVal, jointSpecs)

		
	}
	if (!is.na(fileName))
		dev.off()
	return()
}


#########################################################################################
##	PLOT AAR FRACTIONS
#########################################################################################
plotAARfractionsVage2 <- function (wholeData, aminoList, ageCol, yPrefix, ySuffix, fileName) {

	xData <- wholeData[!is.na(wholeData[, ageCol]),]

	DROP <- "Bush et al."
	c14.specs <- unique(xData$specimen_no)
	for (s in 1:length(c14.specs)) {
		
		c14.spec.data <- xData[(xData[,'specimen_no']==c14.specs[s]),]
		specMethods <- unique(c14.spec.data[,'lab_method'])

		if ((length(specMethods) > 1)) {
			db.na <- is.na(xData[,'lab_method'])
			db.lm <- (xData[,'lab_method']!=DROP)
			db.sp <- (xData[,'specimen_no']!=c14.specs[s])
			xData <- xData[(db.lm | db.sp),]
		}
	}
	
	xData[,ageCol] <-(xData[,ageCol]) + (2013 - 1950)
	
	
	yDataBF <- wholeData[(wholeData[,'analysis_fraction'] == 'B-FAA'),]
	yDataBF  <- useMedianHPLCrunValue(yDataBF, paste(yPrefix, aminoList,ySuffix,sep=''), 'specimen_no')
	jointSpecs <- intersect(xData$specimen_no, yDataBF$specimen_no)

	yDataBT <- wholeData[(wholeData[,'analysis_fraction'] == 'B-THAA'),]
	yDataBT  <- useMedianHPLCrunValue(yDataBT, paste(yPrefix, aminoList,ySuffix,sep=''), 'specimen_no')
	
	yDataAT <- wholeData[(wholeData[,'analysis_fraction'] == 'THAA'),]
	yDataAT  <- useMedianHPLCrunValue(yDataAT, paste(yPrefix, aminoList,ySuffix,sep=''), 'specimen_no')
	
	yLab <- 'D/L'
	xLab <- 'Age (yr)'
	
		
	
	if (!is.vector(jointSpecs))
		return(-1)
	
	xMin <- min(xData[, ageCol], na.rm=TRUE)
	xMax <- max(xData[, ageCol], na.rm=TRUE)
	xMin <- floor((xMin)/100)*100
	xMax <- ceiling((xMax)/100)*100
	xMin <- 5
	xLim <- c(xMin,xMax)
	
	cVal <- paste(yPrefix,aminoList,ySuffix,sep='')

	yMin1 <- min(yDataBF[, cVal], na.rm=TRUE)
	yMax1 <- max(yDataBF[, cVal], na.rm=TRUE)
	yMin2 <- min(yDataBT[, cVal], na.rm=TRUE)
	yMax2 <- max(yDataBT[, cVal], na.rm=TRUE)
	yMin3 <- min(yDataAT[, cVal], na.rm=TRUE)
	yMax3 <- max(yDataAT[, cVal], na.rm=TRUE)
	yMin <- floor(min(yMin1, yMin2, yMin3)*100)/100
	yMax <- ceiling(max(yMax1, yMax2, yMin3)*10)/10
	yLim <- c(yMin,yMax)

	if (!is.na(fileName))
		pdf(paste(PATH.OUT,fileName,sep='/'), width=pageWidthTwo, height=pageWidthTwo)

	par(mfcol=c(ceiling(length(aminoList)/2),2), cex=1, mar=c(4,4,1,1), mgp=c(2,0.7,0), cex.lab=1.25, las=1)
	yLT <- aminoList[seq(1, ceiling(length(aminoList)/2), by=1)]
	xLT <- aminoList[c(ceiling(length(aminoList)/2), length(aminoList))]
	
	for (aa in aminoList) {
	
		DL <- paste(yPrefix,aa,ySuffix,sep='')
		
		ylab <- ''
		if (aa %in% yLT)
			ylab <- yLab
		xlab <- ''
		if (aa %in% xLT)
			xlab <- xLab
		
		if (aa == 'Asp')
			par(mar=c(2,3,1,0))
		if (aa == 'Ala')
			par(mar=c(2,2.5,1,0.5))
		if (aa == 'Glu')
			par(mar=c(3,3,0,0))
		if (aa == 'Phe')
			par(mar=c(3,2.5,0,0.5))
		
		
		plot(1:1, type='n', ylim=yLim,xlim=xLim, ylab=ylab, xlab=xlab)
		text(xMin, yMax*0.9 ,paste(aa,'D/L'), pos=4, cex=1.25)

		yData <- yDataAT

		xyPairs <- c(0,0)

		for (s in jointSpecs) {
			xSpec <- xData[(xData[,'specimen_no']==s),]
			ySpec <- yData[(yData[,'specimen_no']==s),]
			
			xSq <- quantile(xSpec[,DL], qSD)
			ySq <- quantile(ySpec[,DL], qSD)
			
			points(xSq[3], ySq[3], pch=19, cex=0.5)
			segments(xSq[3], ySq[1], xSq[3], ySq[5], lwd=0.5)
			segments(xSq[1], ySq[3], xSq[5], ySq[3], lwd=0.5)
			
			xyPairs <- rbind(xyPairs,c(xSq[3], ySq[3]))
		}

		xyPairs <- xyPairs[-1,]
		
		lm.p <- lm(xyPairs[,2] ~ xyPairs[,1])
		abline(lm.p, lty=2)
		
		text(xMax,0.17, substitute(r^2==k, list(k=round(summary(lm.p)$r.squared,2))), pos=2)
		text(xMax,0.10,paste('slope =',round(summary(lm.p)$coefficients[2,1],2),'\u00b1',round(summary(lm.p)$coefficients[2,2],2)), pos=2)
		text(xMax,0.03,paste('intercept =',round(summary(lm.p)$coefficients[1,1],2),'\u00b1',round(summary(lm.p)$coefficients[1,2],2)), pos=2)	
	}
	
	if (!is.na(fileName))
		dev.off()
	return()
}


#########################################################################################
##	AVERAGE SPECIMEN DATA
#########################################################################################
aarSpecimenAverage <- function (dataOne, iniA3) {

	specs <- unique(dataOne[,'specimen_no'])

	for (s in specs) {
		specRows <- dataOne[(dataOne[,'specimen_no']==s),]
		if (nrow(specRows)>1) {
			for (a in iniA3) {
				specRows[1,paste(a,'DL',sep='_')] <- mean(specRows[,paste(a,'DL',sep='_')],na.rm=TRUE)
			}	
		}
		if (s == specs[1])
			dataSpec <- dataOne[1,]
		else 
			dataSpec <- rbind(dataSpec,specRows[1,])
	}

	return(dataSpec)
}

#########################################################################################
##	DOWNCORE AGE IN D/L
#########################################################################################
plotAARDownCore <- function (dataOne, iniA3, fileName) {

	if (!is.na(fileName))
		pdf(paste(PATH.OUT,fileName,sep='/'), width=pageWidthTwo, height=pageWidthTwo)

	par(mfcol=c(1,length(iniA3)), cex=0.5, mar=c(2,1,2,1), mgp=c(1,0.5,0))
		
	TOPS <- unique(as.factor(areaRaw[,'top']))

	xLab <- expression(DL^{2})

	for (d in iniA3) {
		
		DL <- paste(d,'DL',sep='_')
		
		plot(dataOne[,DL]^2, dataOne[,'top'], ylim=c(170,0), ylab='uncorrected depth (m)', xlab=xLab, main=d, type='n', yaxt='n')
		axis(2,seq(0,180,by=30), seq(0,1.8,by=0.3), las=1)
	
		for (d in TOPS) {
			dataLayer <- dataOne[((dataOne[,'top']==d) & (!is.na(dataOne[,DL]))),]
			if (nrow(dataLayer)>0) {
				a <- quantile(dataLayer[,DL]^2, c(0.025,0.25,0.5,0.75,0.975), rm.na=TRUE)
			#	lines(c(a[1],a[5]),c(d,d))
				bot <- as.numeric(d)-2
				top <- as.numeric(d)+2
				polygon(c(a[2],a[4],a[4],a[2]),c(bot,bot,top,top), col='lightgrey')
			#	points(a[3],d,pch=22, col='blue', bg='blue')
				lines(c(a[3],a[3]),c(top+2,bot-2), lwd=2)
			}
		}
	
		points(dataOne[,DL]^2, dataOne[,'top'], pch=21,col='white',bg='black')
		text(min(dataOne[,DL]^2, na.rm=TRUE),max(dataOne[,'top'], na.rm=TRUE)+10,taxon,font=4, pos=4)
	}
	
	if (!is.na(fileName))
		dev.off()
}


