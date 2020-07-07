## AVERAGE MULTIPLE HPLC RUNS ON THE SAME SPECIMEN
###################################################
mak.calc.Y <- function(dataSet, colx, colu, Y) {
	
	cols <- colnames(dataSet)
	
	sd400 <- (dataSet[,colu] < (Y*1000))
	cv400 <- ((dataSet[,colu]/dataSet[,colx]) < Y)
	use400 <- sd400 | cv400
	
	dataSet <- cbind(dataSet,use400)
	newName <- paste('use',(Y*1000),sep='')
	colnames(dataSet) <- c(cols,newName)

	return(dataSet)
}


## AVERAGE MULTIPLE HPLC RUNS ON THE SAME SPECIMEN
###################################################
averageMultipleHPLC <- function(dataToUse,colR, colL, colD, colFixed) {
		
	specimens <- levels(as.factor(dataToUse[,'specimen_no']))
	columns <- colnames(dataToUse)
	
	dataCon <- matrix(nrow=length(specimens), ncol=length(columns))
	colnames(dataCon) <- columns
	dataCon <- data.frame(dataCon)
	
	for (i in 1:length(specimens)) {
		tmp <- subset(dataToUse, specimen_no==specimens[i])
		if (nrow(tmp) > 1) {
			dataCon[i,] <- tmp[1,]
			for (dl in colR)
				dataCon[i,dl] <- mean(tmp[,dl],na.rm=TRUE)
			for (cL in colL)
				dataCon[i,cL] <- mean(tmp[,cL],na.rm=TRUE)
			for (cD in colD)
				dataCon[i,cD] <- mean(tmp[,cD],na.rm=TRUE)
		} else {
			dataCon[i,] <- tmp
		}
		for (f in colFixed)
			dataCon[i,f] <- as.character(tmp[1,f])
#		dataCon[i,'material'] <- as.character(tmp[1,'material'])
#		dataCon[i,'ual_num'] <- as.character(tmp[1,'ual_num'])
#		dataCon[i,'ansto_num'] <- as.character(tmp[1,'ansto_num'])
	}

	return(dataCon)
}




#
#	SIMPLE CALIBRATION CURVE FOR KNOWN EXPONENT
###################################################
aaCalCurve <- function(j, ageMin, ageMax, data14C, dataAmino, nameAmino, nameSpecimen) {

	if (missing(j))
		j<- 1

	dReg <- lm(data14C ~ I(dataAmino^j))

	ageAmi <- data.frame(dataAmino=sort(dataAmino))
	age.pred <- predict(dReg, int='p', newdata=ageAmi)
	age.conf <- predict(dReg, int='c', newdata=ageAmi)

	if (missing(ageMin)) 
		ageMin <- min(age.pred)

	if (missing(ageMax)) 
		ageMax <- max(age.pred)
	
		xl <- paste(nameAmino,'D/L',j)
		yl <- 'Calibrated Years BP'
		ti <- paste(nameSpecimen,': ',nameAmino,sep='')
		
		plot(1:1,main=ti, ylim=c(ageMin,ageMax), xlim=c(min(dataAmino^j),max(dataAmino^j)), xlab=xl, ylab=yl, pch=19)
		matlines(ageAmi^j,age.pred, lty=c(5,1,1), col=c('black','grey66', 'grey66'))
		matlines(ageAmi^j,age.conf, lty=c(1,1), col=c('grey33', 'grey33'))
		text((min(dataAmino)^j),(ageMax*.9),cex=1.2,pos=4,substitute('D/L'^k, list(k = round(j,2))))
		text((min(dataAmino)^j),(ageMax*.8),cex=1.2,pos=4,substitute(r^2 == k, list(k = round(summary(dReg)$adj.r.squared,3))))
		text((min(dataAmino)^j),(ageMax*.7),cex=1.2,pos=4,paste('max =',round(max(age.pred[,'fit']),0)))
		text((min(dataAmino)^j),(ageMax*.6),cex=1.2,pos=4,paste('min =',round(min(age.pred[,'fit']),0)))
		
		points(ageAmi^j, age.pred[,'fit'], col='blue', pch=19, cex=.8)
		points(dataAmino^j,data14C, col='black', pch=19, cex=.8)

	ageAmi <- data.frame(dataAmino=(dataAmino))
	age.pred <- predict(dReg, int='p', newdata=ageAmi)
	return (age.pred)		
}

#
###################################################
optimAgeSimilData <- function() {
	data14C <- tellina[,'medianAge']
	e14Cold <- tellina[,'X1sOld']
	e14Cyng <- tellina[,'X1sYng']
	dataAA  <- tellina[,p[1]]
	dataGA  <- tellina[,p[2]]
	save(data14C, e14Cold, e14Cyng, dataAA, dataGA, file ='tmp1.Rimage')
}

#
###################################################
optimAgeSimilarity <- function(exps) {

	load('tmp1.Rimage')
	
	expAA <- exps[1]
	expGA <- exps[2]
	
	fitAA <- lm(data14C ~ I(dataAA^expAA), weights=1/(abs(e14Cold-e14Cyng)))
	ageAA <- data.frame(dataAA=sort(dataAA))
	preAA <- predict(fitAA, int='p', newdata=ageAA)

	fitGA <- lm(data14C ~ I(dataGA^expGA), weights=1/(abs(e14Cold-e14Cyng)))
	ageGA <- data.frame(dataAA=sort(dataGA))
	preGA <- predict(fitGA, int='p', newdata=ageGA)

 	fitTT <- lm(preGA[,'fit'] ~ preAA[,'fit'] + 0)
	preTT <- predict(fitGA, int='p')
	
	preTT.var <- preTT[,'upr'] - preTT[,'lwr']
	return(median(preTT.var))
#	return(1-summary(fitTT)$adj.r.squared)
}

#
###################################################
minimzePredInterval <- function(j, ...) {

	load('tmp.Rimage')
	
#	dReg <- lm(data14C[,'X2sMedian'] ~ I(dataAmino^j), weights=1/e14Cspan)

	dReg <- lm(data14C[,'X2sMedian'] ~ I(dataAmino^j) -1, weights=1/e14Cspan)

#	dReg <- lm((data14C[,'X2sMedian']+56) ~ I((dataAmino-min(dataAmino))^j) - 1, weights=1/e14Cspan)

#	ageAmi <- data.frame(dataAmino=sort(dataAmino))
#	age.pred <- predict(dReg, int='p', newdata=ageAmi)
#	age.conf <- predict(dReg, int='c', newdata=ageAmi)

#	pred.var <- age.pred[,'upr'] - age.pred[,'lwr']

#	in weighted regressions this may not be what I want...
#	this may assign too much weight low in the curve...
#	return(median(pred.var))

	return(1 - (summary(dReg))$adj.r.squared)

#	seems to be prone to spectacular failures
#	return(summary(dReg)$coefficients[2])
}

#
###################################################
minimzeRsquared <- function(j, ...) {

	load('tmp.Rimage')
	
#	if (DATEZERO < 2000) {
#		dReg <- lm(data14C[,'X2sMedian'] ~ I(dataAmino^j), weights=(1/e14Cspan)^2)
#	} else {
		dReg <- lm(data14C[,'X2sMedian'] ~ I(dataAmino^j) - 1, weights=(1/e14Cspan)^2)
#	}

	return(1 - (summary(dReg))$adj.r.squared)

}


###################################################
aaCalCurveOptimExponent <- function(j, data14C, weight, dataAmino, DATEZERO) {

	# SET - DEFAULT PARAMETERS
	if (missing(j))
		j<- 1

	# DETERMINE AGE SPAN TO USE FOR WEIGHTING
	e14Cspan <- ((data14C[,'X2sOld']+100) - (data14C[,'X2sYng']+100))

	# DETERMINE MEAN SPAN FOR NON-MODERN SHELLS
	meanSpan <- mean(e14Cspan[(data14C[,'X2sMedian'] > -50)], na.rm=TRUE)
	
	# GIVE EACH MODERN SHELL THE MEAN SPAN FOR NON-MODERN SHELLS
	e14Cspan[(data14C[,'X2sMedian'] < -50)] <- meanSpan
	
	save(data14C, e14Cspan, dataAmino, DATEZERO, file ='tmp.Rimage')
	
	# OPTIMIZE EXPONENT FIT TO DATA
	optPara <- optimize(minimzeRsquared, interval=c(0,6))	
	
	save(data14C, e14Cspan, dataAmino, optPara, file ='tmp.Rimage')

	j <- optPara$minimum

	return(j)
}

#
#	OPTIMIZE REGRESSION LINE
###################################################
aaCalCurveOptim <- function(j, ageMin, ageMax, data14C, weight, dataAmino, nameAmino, nameSpecimen, DATEZERO, LIVING) {

	fontSize <- 1.0
	predPointCol <- 'grey33'
	calPtSize <- 1.1

	if (missing(DATEZERO))
		DATEZERO <- 'XXXX'

	if (missing(LIVING))
		LIVING <- 'XXXX'

	# SET - DEFAULT PARAMETERS
	if (missing(j))
		j<- 1

	if (missing(weight))
		weight <- TRUE

	j <- aaCalCurveOptimExponent(j, data14C, weight, dataAmino, DATEZERO)
	
	load('tmp.Rimage')
	
	##
	## follwing equations needs to match the equation in: 
	##   aaCalCurveOptimExponent() / minimizePredInterval
	##
	
#	if (DATEZERO < 2000) {
#		dReg <- lm(data14C[,'X2sMedian'] ~ I(dataAmino^j), weights=(1/e14Cspan)^2)
#	} else {
		dReg <- lm(data14C[,'X2sMedian'] ~ I(dataAmino^j) - 1, weights=(1/e14Cspan)^2)
#	}
	
	print(summary(dReg)$call)

	ageAmi <- data.frame(dataAmino=sort(dataAmino))
	age.pred <- predict(dReg, int='p', newdata=ageAmi)
	age.conf <- predict(dReg, int='c', newdata=ageAmi)
	
#	if (DATEZERO != 1950) {
#		zerofix <-  DATEZERO - 1950
#		age.pred <- age.pred - zerofix
#	} else {
		zerofix <- 0
#	}

	# SET - PLOT PARAMETERS - AGEMIN & AGEMAX
	if (missing(ageMin)) 
		ageMin <- min(age.pred,na.rm=TRUE)
	if (ageMin > min(data14C[,'X2sYng'],na.rm=TRUE))
		ageMin <- min(data14C[,'X2sYng'],na.rm=TRUE)		
	
	if (missing(ageMax) | is.na(ageMax)) 
		ageMax <- max(age.pred,na.rm=TRUE)
	if (ageMax < max(data14C[,'X2sOld'],na.rm=TRUE))
		ageMax <- max(data14C[,'X2sOld'],na.rm=TRUE)

	xlr <- paste(nameAmino,'D/L')
	xl <- substitute(x ^ k, list(x = xlr, k = round(j,3)))
	yl <- '14C age (years)'
#	yl <- substitute(a^14 + x, list(a = '', x = ylr))
	ti <- paste(nameSpecimen,': ',nameAmino,sep='')

	
	# PLOT - CALIBRATION CURVE
	plot(1:1,main='', ylim=c(ageMin,ageMax), xlim=c(min(dataAmino^j),max(dataAmino^j)), xlab=xl, ylab=yl, pch=19)
#	plot(1:1,main=ti, ylim=c(ageMin,ageMax), xlim=c(min(dataAmino^j),max(dataAmino^j)), xlab=xl, ylab=yl, pch=19)

	# PLOT - ADD CONFIDENCE AND PREDICTION INTERVALS
	matlines(ageAmi^j,age.pred, lty=c(1,1,1), col=c('black','grey66', 'grey66'), lwd=2)
#	matlines(ageAmi^j,age.conf, lty=c(1,1), col=c('grey33', 'grey33'))

	uncert <- age.pred[,'upr'] - age.pred[,'fit']

	# MAKE SUMMARY DATA
	summaryData <- matrix(, nrow=1, ncol=11)
	colnames(summaryData) <- c('dateZero','zeroMod','num','C14','b','m','exp','ar2','unc','min','max')
	rownames(summaryData) <- c(paste(LIVING, weight, nameAmino, nameSpecimen))

	summaryDataN <- matrix(, nrow=1, ncol=4)
	colnames(summaryDataN) <- c('LIVING', 'weight', 'nameAmino', 'nameSpecimen')
	summaryDataN[1,'LIVING'] <- LIVING
	summaryDataN[1,'weight'] <- weight
	summaryDataN[1,'nameAmino'] <- nameAmino
	summaryDataN[1,'nameSpecimen'] <- nameSpecimen

	summaryData[1,'dateZero'] <- DATEZERO
	summaryData[1,'zeroMod'] <- zerofix
	summaryData[1,'C14'] <- nrow(data14C[!is.na(data14C[,'X2sMedian']),])
	summaryData[1,'num'] <- nrow(age.pred)
	summaryData[1,'b'] <- dReg$coefficients['(Intercept)']
	summaryData[1,'m'] <- dReg$coefficients['I(dataAmino^j)']
	summaryData[1,'exp'] <- j
	summaryData[1,'ar2'] <- summary(dReg)$adj.r.squared
	summaryData[1,'unc'] <- round(median(uncert),0)
	summaryData[1,'min'] <- round(min(age.pred[,'fit']),0)
	summaryData[1,'max'] <- round(max(age.pred[,'fit']),0)

	# PLOT - ADD SUMMARY STATS TO PLOT
	text((min(dataAmino)^j),(ageMax*.92),cex=fontSize,pos=4,paste('x =',round(summaryData[1,'exp'],3)))
	text((min(dataAmino)^j),(ageMax*.82),cex=fontSize,pos=4,paste('unc =',summaryData[1,'unc']))
	text((min(dataAmino)^j),(ageMax*.72),cex=fontSize,pos=4,paste('max =',summaryData[1,'max']))
	text((min(dataAmino)^j),(ageMax*.62),cex=fontSize,pos=4,paste('min =',summaryData[1,'min']))
	text((min(dataAmino)^j),(ageMax*.52),cex=fontSize,pos=4,substitute(r^2 == k, list(k = round(summaryData[1,'ar2'],3))))

	text((max(dataAmino)^j),(ageMin),cex=fontSize,pos=2,nameSpecimen, font=3, offset= -0.25)
		
	# PLOT - ADD PREDICTED SAMPLE POINTS
	points(ageAmi^j, age.pred[,'fit'], col=predPointCol, pch=21, cex=.8)

	# PLOT - ADD CARBON DATED SHELLS W/ ERROR BARS
	points(dataAmino^j,data14C[,'X2sMedian'], col='black', pch=22, cex=calPtSize, lwd=1)
	for (i in 1:length(data14C[,'X2sYng']))
		lines(c(dataAmino[i]^j,dataAmino[i]^j),c(data14C[i,'X2sOld'], data14C[i,'X2sYng']), lty=1, lwd=1)

	# WRITE SUMMARY DATA
	summaryData <- cbind(summaryDataN,summaryData)
	fileConnection <- paste(FILEOUT, taxon,'-calibSummary.csv',sep='')
	write.table(summaryData, fileConnection, append = TRUE, col.names=TRUE)

	# RETURN CALIBRATED AGES (I.E. PREDICTED AGES)
	ageAmi <- data.frame(dataAmino=(dataAmino))
	age.pred <- predict(dReg, int='p', newdata=ageAmi)
	return (age.pred)		
}


#
#	SIMPLE CALIBRATION CURVE WITH ERROR
###################################################
aaCalCurveError <- function(j, ageMin, ageMax, data14C, e14Cold, e14Cyng, dataAmino, nameAmino, nameSpecimen) {

	if (missing(j))
		j<- 1

	dReg <- lm(data14C ~ I(dataAmino^j), weights=1/(abs(e14Cold-e14Cyng)))

	ageAmi <- data.frame(dataAmino=sort(dataAmino))
	age.pred <- predict(dReg, int='p', newdata=ageAmi)
	age.conf <- predict(dReg, int='c', newdata=ageAmi)

	if (missing(ageMin)) 
		ageMin <- min(age.pred)

	if (missing(ageMax)) 
		ageMax <- max(age.pred)
	
		xl <- paste(nameAmino,'D/L',round(j,3))
		yl <- 'Calibrated Years BP'
		ti <- paste(nameSpecimen,': ',nameAmino,sep='')
		
		plot(1:1,main=ti, ylim=c(ageMin,ageMax), xlim=c(min(dataAmino^j),max(dataAmino^j)), xlab=xl, ylab=yl, pch=19)
		matlines(ageAmi^j,age.pred, lty=c(5,1,1), col=c('black','grey66', 'grey66'))
		matlines(ageAmi^j,age.conf, lty=c(1,1), col=c('grey33', 'grey33'))
		text((min(dataAmino)^j),(ageMax*.9),cex=1.2,pos=4,substitute('D/L' == k, list(k = round(j,3))))
		text((min(dataAmino)^j),(ageMax*.8),cex=1.2,pos=4,substitute(r^2 == k, list(k = round((summary(dReg))$adj.r.squared,3))))
		text((min(dataAmino)^j),(ageMax*.7),cex=1.2,pos=4,paste('unc =',round(median(age.pred[,'upr']-age.pred[,'lwr']),0)))
		text((min(dataAmino)^j),(ageMax*.6),cex=1.2,pos=4,paste('max =',round(max(age.pred[,'fit']),0)))
		text((min(dataAmino)^j),(ageMax*.5),cex=1.2,pos=4,paste('min =',round(min(age.pred[,'fit']),0)))
		
		points(ageAmi^j, age.pred[,'fit'], col='blue', pch=19, cex=.8)
		points(dataAmino^j,data14C, col='black', pch=19, cex=.8)
		for (i in 1:length(e14Cold))
			lines(c(dataAmino[i]^j,dataAmino[i]^j),c(e14Cold[i], e14Cyng[i]),cex=0.1, lty=1)

	ageAmi <- data.frame(dataAmino=(dataAmino))
	age.pred <- predict(dReg, int='p', newdata=ageAmi)
	return (age.pred)		
}

#
#	SIMPLE AGE HISTOGRAM
###################################################
aaCalHist <- function(dAge, ageMin, ageMax, breaks, nameAmino, nameSpecimen) {

	if (missing(ageMin)) 
		ageMin <- min(dAge[,2])

	if (missing(ageMax)) 
		ageMax <- max(dAge[,2])

	xl <- 'calibrated YBP'
	yl <- 'specimens'
	ti <- paste(nameSpecimen,': ',nameAmino,sep='')

	hist(dAge[,2], xlim=c(ageMin,ageMax), breaks=breaks, freq=TRUE, xlab=xl, ylab=yl, main=ti, col='grey')
}
