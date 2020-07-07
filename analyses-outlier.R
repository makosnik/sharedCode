library(outliers)

mak.wrapper.wt.calcs <- function(dataSet,aminos) {

	big <- 10000
	
	ageCols <- vector(mode='logical', length(aminos))
	for (amino in 1:length(aminos))
		ageCols[amino] <- paste(aminos[amino],'AgeFit',sep='_')
	
	uncCols <- vector(mode='logical', length(aminos))
	for (amino in 1:length(aminos))
		uncCols[amino] <- paste(aminos[amino],'AgeUnc',sep='_')

	ageUncert <- matrix(nrow=nrow(dataSet),ncol=length(aminos))
	colnames(ageUncert)<-uncCols

	medAll <- matrix(nrow=nrow(dataSet),ncol=4)
	colnames(medAll) <- c('medLwr','medMed','medUpr','medUnc')

	g10.p <- vector(mode='logical',nrow(dataSet))
	g10.h <- vector(mode='logical',nrow(dataSet))
	g10.wm <- vector(mode='logical',nrow(dataSet))
	g10.wu <- vector(mode='logical',nrow(dataSet))

	g20.p <- vector(mode='logical',nrow(dataSet))
	g20.h <- vector(mode='logical',nrow(dataSet))

	g11.p <- vector(mode='logical',nrow(dataSet))
	g11.h <- vector(mode='logical',nrow(dataSet))

	gMinus <- vector(mode='logical',nrow(dataSet))

	wmAll <- vector(mode='logical',nrow(dataSet))
	wuAll <- vector(mode='logical',nrow(dataSet))
	sdAll <- vector(mode='logical',nrow(dataSet))
	emAll <- vector(mode='logical',nrow(dataSet))

	wmM1 <- vector(mode='logical',nrow(dataSet))
	wuM1 <- vector(mode='logical',nrow(dataSet))
	sdM1 <- vector(mode='logical',nrow(dataSet))
	emM1 <- vector(mode='logical',nrow(dataSet))
	minus1 <- vector(mode='logical',nrow(dataSet))

	wmM2 <- vector(mode='logical',nrow(dataSet))
	wuM2 <- vector(mode='logical',nrow(dataSet))
	sdM2 <- vector(mode='logical',nrow(dataSet))
	emM2 <- vector(mode='logical',nrow(dataSet))
	minus2 <- vector(mode='logical',nrow(dataSet))

	wmB1 <- vector(mode='logical',nrow(dataSet))
	wuB1 <- vector(mode='logical',nrow(dataSet))
	sdB1 <- vector(mode='logical',nrow(dataSet))
	emB1 <- vector(mode='logical',nrow(dataSet))
	minusB <- vector(mode='logical',nrow(dataSet))


	for (amino in 1:length(aminos))
		ageUncert[,paste(aminos[amino],'AgeUnc',sep='_')] <- (dataSet[,paste(aminos[amino],'AgeUpr',sep='_')]+big) - (dataSet[,paste(aminos[amino],'AgeLwr',sep='_')]+big)
	
	## CANNOT ALLOW AN UNCERTAINTY OF 0
	ageUncert[(ageUncert[,]==0)] <- 1

	for (i in 1:nrow(dataSet)) {
		
		ageCol4 <- ageCol3 <- ageCol2 <- ageCols
		uncCol4 <- uncCol3 <- uncCol2 <- uncCols

		## NON-PARAMETRIC ESTIMATE
#		medAll[i,c('medLwr','medMed','medUpr')] <- unlist(quantile(dataSet[i,ageCols],c(0.159,0.5,0.841)))
		medAll[i,c('medLwr','medMed','medUpr')] <- unlist(quantile(dataSet[i,ageCols],c(0.25,0.5,0.75)))
		medAll[i,'medUnc'] <- (medAll[i,'medUpr']+1000)-(medAll[i,'medLwr']+1000)

		## ALL SPECIMENS - WEIGHTED ESTIMATE
		wmAll[i] <- round(mak.wt.mean(unlist(dataSet[i,ageCols]),unlist(1/ageUncert[i,uncCols]^2)),0)
		sdAll[i] <- round(mak.wt.sdev(unlist(dataSet[i,ageCols]),unlist(1/ageUncert[i, uncCols]^2)),0)
		emAll[i] <- round(mak.wt.errorOfMean(unlist(1/ageUncert[i, uncCols]^2)),0)
		wuAll[i] <- sdAll[i]
		if (emAll[i] > wuAll[i])
			wuAll[i] <- emAll[i]

		mdAll <- unlist(abs(dataSet[i,ageCols] - wmAll))
				
		if (length(ageCols) > 5) {
			
			for (a in 1:length(ageCols))
				if (mdAll[a] == max(mdAll)) {
					minus1[i] <- substr(ageCols[a],1,3)
					ageCol3 <- setdiff(ageCols,c(paste(minus1[i],'AgeFit',sep='_')))
					uncCol3 <- setdiff(uncCols,c(paste(minus1[i],'AgeUnc',sep='_')))
			}
			
			wmM1[i] <- round(mak.wt.mean(unlist(dataSet[i,ageCol3]),unlist(1/ageUncert[i,uncCol3]^2)),0)
			sdM1[i] <- round(mak.wt.sdev(unlist(dataSet[i,ageCol3]),unlist(1/ageUncert[i, uncCol3]^2)),0)
			emM1[i] <- round(mak.wt.errorOfMean(unlist(1/ageUncert[i, uncCol3]^2)),0)
			wuM1[i] <- sdM1[i]
			if (emM1[i] > wuM1[i])
				wuM1[i] <- emM1[i]

			mdM1 <- unlist(abs(dataSet[i,ageCol3]-wmM1))

			for (a in 1:length(ageCol3))
				if (mdM1[a] == max(mdM1)) {
					minus2[i] <- substr(ageCol3[a],1,3)
					ageCol4 <- setdiff(ageCol3,c(paste(minus2[i],'AgeFit',sep='_')))
					uncCol4 <- setdiff(uncCol3,c(paste(minus2[i],'AgeUnc',sep='_')))
			}
			
			wmM2[i] <- round(mak.wt.mean(unlist(dataSet[i,ageCol4]),unlist(1/ageUncert[i,uncCol4]^2)),0)
			sdM2[i] <- round(mak.wt.sdev(unlist(dataSet[i,ageCol4]),unlist(1/ageUncert[i, uncCol4]^2)),0)
			emM2[i] <- round(mak.wt.errorOfMean(unlist(1/ageUncert[i, uncCol4]^2)),0)				
			wuM2[i] <- sdM2[i]
			if (emM2[i] > wuM2[i])
				wuM2[i] <- emM2[i]
			
#			g20 <- grubbs.test(unlist(dataSet[i,ageCols]), type=20, two.sided=TRUE)
#			g20.p[i] <- g20$p.value			
#			g20.h[i] <- g20$alternative			
			g10 <- grubbs.test(unlist(dataSet[i,ageCols]), type=10, two.sided=TRUE)
#			g11 <- grubbs.test(unlist(dataSet[i,ageCols]), type=11, two.sided=TRUE)
#			g11.p[i] <- g11$p.value			
#			g11.h[i] <- g11$alternative			

			g10.p[i] <- g10$p.value			
			g10.h[i] <- g10$alternative			
			
			pCrit <- (0.05/nrow(dataSet))
						
			if (g10$p.value < pCrit) {
				killName <- substr(unlist(dimnames(as.array(g10$statistic[1]))),3,5)
				gMinus[i] <- killName
				dd <- mak.wt.calcs(dataSet[i,],setdiff(aminos,c(killName)))
				g10.wm[i] <- dd[1,'wmAll']
				g10.wu[i] <- dd[1,'wuAll']
			}
		}
	}
	
#	dataSet <- cbind(dataSet,g10.p,g10.h,g20.p,g20.h,g11.p,g11.h, gMinus)
	dataSet <- cbind(dataSet,wmAll,wuAll,sdAll,emAll,wmM1,wuM1,sdM1,emM1,minus1,wmM2,wuM2,sdM2,emM2,minus2,medAll,g10.p,g10.h,gMinus,g10.wm,g10.wu)
	
	return(dataSet)
}


mak.wrapper.plot.ages <- function(dataToUse, ymin, ymax) {

	iniA <- c('Asp','Glu','Ala','Val','Phe','Leu','A.I')

	ageUpr <- vector(mode='logical', length(iniA))
	ageLwr <- vector(mode='logical', length(iniA))
	ageFit <- vector(mode='logical', length(iniA))
	
	for (amino in 1:length(iniA)) {
		ageLwr[amino] <- paste(iniA[amino],'AgeLwr',sep='_')
		ageUpr[amino] <- paste(iniA[amino],'AgeUpr',sep='_')
		ageFit[amino] <- paste(iniA[amino],'AgeFit',sep='_')
				
	}

alab <- 0.75

i <- 1
for (i in 1:nrow(dataToUse)) {
	
	if (missing(ymax))
		ymax <- ceiling((max(dataToUse[i,c(ageUpr,'X2sOld')], na.rm=TRUE)*1.03)/100)*100
	if (missing(ymin))
		ymin <- floor((min(dataToUse[i,c(ageLwr,'X2sYng')], na.rm=TRUE)*0.97)/100)*100


	plot(1:1,type='n',ylim=c(ymin,ymax),xlim=c(0,13), xaxt='n', ylab='Age (a BP)', xlab='',xaxs='i',yaxs='i')
	title(paste(dataToUse[i,'material'],', specimen #',dataToUse[i,'specimen_no'],sep=''),font.main=1)
	
	for (a in 1:length(iniA)) {
		points(a,dataToUse[i,paste(iniA[a],'AgeFit',sep='_')], pch=19)
		lines(c(a,a), c(dataToUse[i, paste(iniA[a],'AgeLwr',sep='_')],dataToUse[i, paste(iniA[a],'AgeUpr',sep='_')]))
		text(a+0.3,dataToUse[i,paste(iniA[a],'AgeLwr',sep='_')], iniA[a], cex= alab)
	}
	
	a <- a+1
	points(a,dataToUse[i,'X2sMedian'], pch=19)
	lines(c(a,a), c(dataToUse[i,'X2sYng'],dataToUse[i,'X2sOld']))
	text(a+0.3, dataToUse[i,'X2sYng'], '14C', cex=alab, pos=1)
	
	a <- a+1
	lines(c(a,a),c(ymin,ymax))
#	axis(1,c(4,10.5),c('individual ages','combined ages'))
	
#	a <- a+1
#	points(a,dataToUse[i,'medMed'], pch=19)
#	lines(c(a,a), c(dataToUse[i,'medLwr'],dataToUse[i,'medUpr']))
#	text(a+0.3,dataToUse[i,'medLwr'], 'NP', cex=alab, pos=1)

#	a <- a+1
#	points(a,dataToUse[i,'umAll'], pch=19)
#	lines(c(a,a), c(dataToUse[i,'umAll']+dataToUse[i,'udAll'],dataToUse[i,'umAll']-dataToUse[i,'udAll']))
#	text(a+0.3,(dataToUse[i,'umAll']-dataToUse[i,'udAll']), 'mean', cex=alab, pos=1)

	a <- a+1
	points(a,dataToUse[i,'wmAll'], pch=19)
	lines(c(a,a), c(dataToUse[i,'wmAll']+ dataToUse[i,'wuAll'], dataToUse[i,'wmAll']-dataToUse[i,'wuAll']))
	text(a+0.3,(dataToUse[i,'wmAll']-dataToUse[i,'wuAll']), 'wmean', cex=alab, pos=1)

	a <- a+1
	dd <- mak.wt.calcs(dataToUse[i,],c('Asp','Glu'))
	points(a,dd[i,'wmAll'], pch=19)
	lines(c(a,a), c(dd[i,'wmAll']+dd[i,'wuAll'],dd[i,'wmAll']-dd[i,'wuAll']))
	text(a+0.3,(dd[i,'wmAll']-dd[i,'wuAll']), 'Asp&Glu', cex=alab, pos=1)

	if (dataToUse[i,'material']=='Liloa') {
		a <- a+1
		dd <- mak.wt.calcs(dataToUse[i,],c('Asp','Ala','Val','Phe','Leu','A.I'))
		points(a,dd[i,'wmAll'], pch=19)
		lines(c(a,a), c(dd[i,'wmAll']+dd[i,'wuAll'],dd[i,'wmAll']-dd[i,'wuAll']))
		text(a+0.3,(dd[i,'wmAll']-dd[i,'wuAll']), '-Glu', cex=alab, pos=1)
	}
	

#	a <- a+1
#	points(a,dataToUse[i,'wmM1'], pch=19)
#	lines(c(a,a), c(dataToUse[i,'wmM1']+ dataToUse[i,'wuM1'], dataToUse[i,'wmM1']-dataToUse[i,'wuM1']))
#	text(a+0.3,dataToUse[i,'wmM1']-dataToUse[i,'wuM1'], paste('-',dataToUse[i,'minus1']), cex=alab, pos=1)

#	a <- a+1
#	points(a,dataToUse[i,'wmM2'], pch=19)
#	lines(c(a,a), c(dataToUse[i,'wmM2']+ dataToUse[i,'wuM2'], dataToUse[i,'wmM2']-dataToUse[i,'wuM2']))
#	text(a+0.3,dataToUse[i,'wmM2']-dataToUse[i,'wuM2'], paste('-',dataToUse[i,'minus2']), cex=alab, pos=1)


	}
}



mak.wrapper.varationPlot <- function (dataCalib,fileName,xCol,yCol,TAXA,maxSD,maxCV) {
	
	print(paste("writing eps file: ", fileName, sep=''))
	
	postscript(fileName)
	par(mar=c(4,4,2,2), las=0, mgp=c(2,1,0))
	plotSize <- 7
	rows <- 3
	cols <- 2
	nplots <- rows*cols
	layout(matrix(seq(nplots), rows, cols, byrow=TRUE), width=rep(lcm(plotSize),cols), height=rep(lcm(plotSize),rows), respect=TRUE)
	plotLetter <- 1
	
#	dataCalib[,xCol] <- dataCalib[,xCol] + 56
	
	if (missing(maxSD))
		maxSD <- max(dataCalib[,yCol])
	if (missing(maxCV))
		maxCV <- max(abs(dataCalib[,yCol]/dataCalib[,xCol]),na.rm=TRUE)
	xmax <- max(dataCalib[,xCol])
	
	for (taxon in TAXA) {
		datesToUse <- dataCalib[(dataCalib[,'material']==taxon),]
	
		plot(datesToUse[,xCol], datesToUse[,yCol], main='', ylab='uncertainty of the mean (a)',xlab='mean inferred age (a BP)', pch=21, cex=1, col='white', bg='black', xlim=c(0,xmax),ylim=c(0,maxSD))

		title (main=paste(LETTERS[plotLetter],'.',sep=''),adj=0)
		title (main=taxon, adj=0.5, font.main=3)
		plotLetter <- plotLetter+1
	
		# Y lines
		for (y in c(0.05,0.1,0.2,0.4)) {
			x <- seq((1000),(xmax*1.3),length.out=2)
			lines(x,x*y, lty=3)
			z <- c(-100,1000)
			lines(z,rep((1000*y),2),lty=3)
		}
		lines(c(xmax*1.3,-100),c(0,0),lty=3)
		
		plot(datesToUse[,xCol],abs(datesToUse[,yCol]/datesToUse[,xCol])+0.01, main='', ylab='coefficient of variation',xlab='mean inferred age (a BP)', pch=21, cex=1, col='white', bg='black', xlim=c(-100,xmax),ylim=c(0.01,100), log='y')

		title (main=paste(LETTERS[plotLetter],'.',sep=''),adj=0)
		title (main=taxon, adj=0.5, font.main=3)
		plotLetter <- plotLetter+1
	
		# Y lines
		for (y in c(0.05,0.1,0.2,0.4)) {
			x <- seq((1000),xmax,length.out=2)
			lines(x,rep(y,2),lty=3)
			z <- seq(10,1000,by=10)
			lines(z,(1000*y)/z, lty=3)
		}
	
	}
	
	dev.off()

}
