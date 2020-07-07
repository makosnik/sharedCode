source('../sharedCode/histogram.R')

# LOAD & CLEAN DATA
###################################################
getCleanAAR <- function(dataRaw,cutQ,taxon,colAS,namAA,fileName) {

	# SUBSET DATA
	#####################################
	if (taxon == 'All')
		dates <- dataRaw
	else 
		dates <- subset(dataRaw, material==taxon)
	
	# OUTPUT REJECTED DATA
	#####################################
	of <- file(paste("./output/variance-",taxon,'-',fileName,'-',cutQ,"-rejectList",".txt",sep=''),'wt')
	sink(of)
	
	
		# CLEAN DATA - REMOVE 0.05 MOST VARIABLE DATA
		#####################################
	
		i<-1
		q1 <- quantile(dates[,colAS[i]],c(cutQ))
		reject1 <- dates[dates[,colAS[i]]>=q1,'specimen_no']
		print(paste(taxon," : rejected shells : ",namAA[i],sep=''))
		print(reject1)
	
		if (length(colAS)>=2) {
			i<-2
			q2 <- quantile(dates[,colAS[i]],c(cutQ))
			reject2 <- dates[dates[,colAS[i]]>=q2,'specimen_no']
			print(paste(taxon," : rejected shells : ",namAA[i],sep=''))
			print(reject2)
		}
		
		if (length(colAS)>=3) {
			i<-3
			q3 <- quantile(dates[,colAS[i]],c(cutQ))
			reject3 <- dates[dates[,colAS[i]]>=q3,'specimen_no']
			print(paste(taxon," : rejected shells : ",namAA[i],sep=''))
			print(reject3)
		}
		
		i<-1
		dates <- dates[dates[,colAS[i]]<q1,]
		if (length(colAS)>=2) {
			i<-2
			dates <- dates[dates[,colAS[i]]<q2,]
		}
		if (length(colAS)>=3) {
			i<-3
			dates <- dates[dates[,colAS[i]]<q3,]
		}
		
	sink()
	
	return(dates)
}

# PLOT BIG VARIANCE PLOT
###################################################
variancePlotsAARbyTaxon <- function(dates, colAA, namAA, ci, colAS, varMeasure) {

	# CONSTANTS
	###################################################
	TAXA <- levels(as.factor(dates[,'material']))

	was3 <- 2
	nplots <- was3*length(TAXA)
	layout(matrix(seq(nplots), length(TAXA), was3, byrow=TRUE), width=rep(2,3), height=rep(2,length(TAXA)), respect=TRUE)

	plotLett <- 0
	
	for (taxon in TAXA) {

		tdates <- subset(dates, material==taxon)

		for (i in 1:length(colAA)) {
		
			# HISTOGRAM - D/L RATIO
#			histCI(tdates[,colAA[i]], xlab=paste(namAA[i],' D/L'), ylab='specimens', title=paste(taxon,' : ',namAA[i],' : D/L'))
			
			# HISTOGRAM - STANDARD DEVIATION ON MEASUREMENT
#			histCI(tdates[,colAS[i]], ci, xlab=paste(namAA[i],varMeasure), ylab='specimens', title=paste(taxon,namAA[i],varMeasure,sep=': '))
		
			# PLOT D/L v. SD
			plotLett <- plotLett+1
#			ti <- paste(LETTERS[plotLett],'. ',taxon,' ',varMeasure,' v. D/L', sep='')
			ti <- paste(LETTERS[plotLett],'. ',taxon, sep='')
			ti <- ''
			plotDLSD (tdates[,c(colAA[i],colAS[i],'specimen_no')], xlab=namAA[i],ylab=paste(namAA[i],varMeasure), main=ti, varMeasure, ci)
			title (main=paste(LETTERS[plotLett],'.',sep=''),adj=0)
			title (main=taxon, font.main=4)
		}
	}
}


# PLOT BIG VARIANCE PLOT
###################################################
variancePlotsAAR <- function(dates,taxon, colAA, namAA, ci, colAS, varMeasure) {

	# SUBSET DATA
	#####################################
	if (taxon != 'All')
		dates <- subset(dates, material==taxon)

	plotSize <- 7
	rows <- 3
	cols <- length(colAA)
	nplots <- rows*cols
	layout(matrix(seq(nplots), rows, cols, byrow=FALSE), width=rep(lcm(plotSize),cols), height=rep(lcm(plotSize),rows), respect=TRUE)
		
	for (i in 1:length(colAA)) {
		
		# HISTOGRAM - D/L RATIO
		L <- i
		title <- paste(LETTERS[L],'. ',namAA[i],' D/L',sep='')
		histCI(dates[,colAA[i]], xlab=' D/L', ylab='number of specimens', title=title)
		
		# HISTOGRAM - STANDARD DEVIATION ON MEASUREMENT
		L <- L+cols
		title <- paste(LETTERS[L],'. ',namAA[i],' ',varMeasure,sep='')
		histCI(dates[,colAS[i]], ci, xlab=varMeasure, ylab='number of specimens', title=title)
	
		# SD - CUMLATIVE DIST
#		plotRank (dates[,colAS[i]], ylab=paste(namAA[i],varMeasure), xlab='shell rank', main=paste(varMeasure,'by rank'), ci)
			
		# PLOT D/L v. SD
		L <- L+cols
		title <- paste(LETTERS[L],'. ',namAA[i],' ',varMeasure,' v. D/L',sep='')
		plotDLSD (dates[,c(colAA[i],colAS[i],'specimen_no')], xlab=paste(namAA[i],' D/L'),ylab=paste(namAA[i],varMeasure), main=title, varMeasure, ci)
	}
}

plotDLSD <- function (pData,xlab,ylab,main,varMeasure,ci,specs) {

		plot (pData[,1],pData[,2], xaxs='i', yaxs='i', xlim=c(0,max(pData[,1], na.rm=TRUE)*1.1), ylim=c(0,max(pData[,2], na.rm=TRUE)*1.1), xlab=xlab, ylab=ylab, type='n', main=main)
		points(pData[,1],pData[,2], pch=21, col='white', bg='black', cex=0.75, lwd=0.25)
		
		# MASK IMPOSSIBLE VALUES
		if (varMeasure=='SD') {
			topOut <- max(pData[,2], rm.na=TRUE)
			topOut <- 0.6
		
			polygon(c(0,topOut,0),c(0,topOut,topOut), col='grey', lty=0)
			lines(c(0,topOut),c(0,topOut*0.50), col='grey')
			lines(c(0,topOut),c(0,topOut*0.25), col='grey')
			
			print(paste('DLSD - grey plotted',topOut,main))
		}
		
		if (!is.logical(ci)) {
			# QUANTILES OF STANDARD DEVIATION ON MEASUREMENT
			hq <- quantile(pData[,2],ci, na.rm=TRUE)
	
			# PLOT SD QUANTILE LINES 
			for (c in 1:length(ci))
				lines(c(0,length(pData[,2])),c(hq[c],hq[c]),lwd=1, lty=2, col='black')

			# GET OFFSET FOR TEXT
			xmin <- min(pData[,2], na.rm=TRUE)
			xmax <- max(pData[,2], na.rm=TRUE)
			xoff <- (xmax - xmin)*0.05

			textPos <- 3
			# PLOT SPECIMEN NUMBERS FOR THE REALLY BAD SPECIMENS 
			for (p in 1:length(pData[,2])) {
#				if (pData[p,3] %in%c(5624,5608,5289,5268,5480,5609,4722,5681)) {
#				if (pData[p,3] %in%c(5624,5608,5289,5268,5377,4824,4722,5681,5609,4732,4761)) {
				if (pData[p,3] %in%c(5624,5608,5289,5268,5377,4824,5335,5681,5609,4772,5455,5676)) {
					text(pData[p,1],pData[p,2],pData[p,3], pos=textPos, offset=0.2, cex=0.8)
					points(pData[p,1],pData[p,2], pch=21, col='grey', bg='grey', lwd=0.25, cex=0.8)
				}
			}

			# PLOT SPECIMEN NUMBERS FOR THE REALLY BAD SPECIMENS 
#			for (p in 1:length(pData[,2])) {
#				if (!is.na(pData[p,2]))
#				if (pData[p,2]>hq[1]) {
#
#					textPos <- 3
#					if (pData[p,3] %in%c(5428,4780,4788,5450))
#						textPos <- 1
#					if (pData[p,3] %in%c(4808,5332,4722,4715,5368,5450,4819,4772))
#						textPos <- 2
#					if (pData[p,3] %in%c())
#						textPos <- 3
#					if (pData[p,3] %in%c(5403,5447,5353,5315,5616,4738,4749,4824,4732,4729,4775,5367,4886,4736))
#						textPos <- 4
#					
#					text(pData[p,1],pData[p,2],pData[p,3], pos=textPos, offset=0.2, cex=0.8)
#				}
#			}
		}		
}

plotRank <- function (rData, xlab, ylab, main, ci) {

	plot (sort(rData), ylab=ylab, xlab=xlab, pch=19, cex=0.5, main=main)

	if (!is.logical(ci)) {
		# QUANTILES OF STANDARD DEVIATION ON MEASUREMENT
		hq <- quantile(rData,ci)
	
		# PLOT SD QUANTILE LINES 
		for (c in 1:length(ci)) {
			lines(c(length(rData)*.1,length(rData)),c(hq[c],hq[c]),lwd=1, lty=2, col='black')
			text((length(rData)*.05),hq[c],ci[c], cex=0.8)
		}
	}
}


aVaPlotsAAR <- function(dates, xNum, yNum, taxon, colAA, namAA, expAA) {

	ylabel <- paste(namAA[yNum],'^',expAA[yNum])
	xlabel <- paste(namAA[xNum],'^',expAA[xNum])
	
	regressPlotOpim(dates[,colAA[xNum]]^expAA[xNum],dates[,colAA[yNum]]^expAA[yNum], yl=ylabel, xl=xlabel, ti=taxon, pl=dates[,'specimen_no'], pun=FALSE, pr2=TRUE)

#	CARBON DATING STUFF

#	axis(side=4,at=tmpY,labels=as.character(tmpD))
	
#	for (i in 1:length(tmpY)) {
#		lines(c(tmpX[i]-1,tmpX[i]),c(tmpY[i],tmpY[i]),col='blue')
#			text(tmpX[i],tmpY[i],tmpD[i],col='blue')
#	}

}