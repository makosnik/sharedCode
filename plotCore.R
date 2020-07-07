#
#	dRef needs to be a matrix with columns X, Y (and total for proportion)
#
plotCoreError<-function(dRef, title, y_lab, y_top, y_bot, x_lab, x_old, x_yng, proportion, logBase, cexMain) {

	# TRAP / SET MISSING VALUES
	if (missing(proportion)) {
		proportion <- FALSE
		if (missing(x_lab))
			x_lab <- '(g)'
		}

	if (proportion) {
		dRef[,2] <- dRef[,2]/dRef[,3]
		if (missing(x_lab))
			x_lab <- 'proportion of layer'
	}

	if (!missing(logBase)) {
		dRef[,2] <- log(logBase,dRef[,2])
		
		if (missing(x_lab))
			x_lab <- paste('log',logBase)
	}
	
	if (missing(y_lab))
		y_lab <- 'core depth (cm)'
	
	if (missing(y_top)) 
		y_top <- min(dRef[,1])

	if (missing(y_bot)) 
		y_bot <- max(dRef[,1])

	if (missing(x_lab))
		x_lab <- ''

	if (missing(x_old))
		x_old <- 0

	if (missing(x_yng))
		x_yng <- 1

	if (missing(cexMain)) {
		cexMain <- 1.6
	}

	# SET CONSTANTS / PREFERENCES
	par(las=1, xaxs="i", yaxs="i", cex.axis=cexMain, cex.main=cexMain, cex.lab=cexMain)	

	tRef <- levels(as.factor(dRef[,'top']))
		
	plot ((dRef[,2]),dRef[,'top'], ylim=c(y_bot,y_top), ylab=y_lab, xlim=c(x_yng, x_old), xlab=x_lab, main=title, cex=0.5)

	for (t in 1:length(tRef)) {
		if (length(dRef[(dRef[,'top']==as.numeric(tRef[t])),2])>=3) {
			lm <- mean((dRef[(dRef[,'top']==as.numeric(tRef[t])),2]), na.rm=T)
			points(lm,tRef[t],cex=1.5, pch=5)
			ls <- sd ((dRef[(dRef[,'top']==as.numeric(tRef[t])),2]), na.rm=T)
			points(lm-ls,tRef[t],cex=1, pch='I')			
			points(lm+ls,tRef[t],cex=1, pch='I')
			lines (c(lm+ls,lm-ls),c(tRef[t],tRef[t]))
		}
	}
}

plotCorePolygonCI<-function(dRef, title, y_lab, y_top, y_bot, x_lab, x_old, x_yng, proportion, division, cis, colours, cexMain) {

	# TRAP / SET MISSING VALUES
	if (missing(y_lab))
		y_lab <- 'core depth (cm)'
	
	if (missing(y_top)) 
		y_top <- min(dRef[,1])

	if (missing(y_bot)) 
		y_bot <- max(dRef[,1])

	if (missing(x_lab))
		x_lab <- ''

	if (missing(x_old))
		x_old <- 0

#	if (missing(x_yng))
#		x_yng <- 1

	if (missing(proportion))
		proportion <- TRUE

	if (missing(division))
		division <- -1

	if (missing(cis))
		cis <- c(0.05,0.25,0.5,0.75,0.95)

	if (missing(colours))
		colours <- c('green','blue','red')

	if (missing(cexMain))
		cexMain <-1.6

	# SET CONSTANTS / PREFERENCES
	par(las=1, xaxs="i", yaxs="i", cex.axis=cexMain, cex.main=cexMain, cex.lab=cexMain)	

	tRef <- levels(as.factor(dRef[,'top']))

	if (proportion) {
		if (x_lab=='x')
			x_lab <- 'proportion of layer'
		dRef[,2] <- dRef[,2]/dRef[,3]
	}

	if (missing(x_yng))
		x_yng <- max(dRef[,2], na.rm=T)*1.1

	# MEDIAN & CIs
	if (division > 0)
		cisTop <- quantile(dRef[(dRef[,'top']<=division),2],cis,rm.na=TRUE)
	cisBot <- quantile(dRef[(dRef[,'top']>division),2],cis,rm.na=TRUE)

	# PLOT AXES
	plot (dRef[,2],dRef[,'top'], ylim=c(y_bot,y_top), ylab=y_lab, xlim=c(x_old,x_yng), xlab=x_lab, main=paste(title), type='n')
	
	# PLOT TOP LAYER - CI & MEDIAN
	if (division > 0) {
		for (i in 1:((length(cisTop)-1)/2)) {
			j <- length(cisTop) +1 - i
			polygon(c(cisTop[i],cisTop[j],cisTop[j],cisTop[i]),c(y_top,y_top,division,division),lty=0, col=colours[i])
		}
		lines(c(cisTop['50%'],cisTop['50%']),c(y_top,division),lwd=3, lty=1, col=colours[i+1])
	}

	# PLOT BOTTOM LAYER - CI & MEDIAN
	for (i in 1:((length(cisBot)-1)/2)) {
		j <- length(cisBot) +1 - i
		polygon(c(cisBot[i],cisBot[j],cisBot[j],cisBot[i]),c(division,division,maxDepth,maxDepth),lty=0, col=colours[i])
	}
	lines(c(cisBot['50%'],cisBot['50%']),c(division,maxDepth),lwd=3, lty=1, col=colours[i+1])
	
	# PLOT POINTS + ERROR BARS
	points (dRef[,2],dRef[,'top'], cex=0.5)
	for (t in 1:length(tRef)) {
		fmin <- min(dRef[(dRef[,'top']==as.numeric(tRef[t])),2], na.rm=T)
		fmax <- max(dRef[(dRef[,'top']==tRef[t]),2], na.rm=T)
		if (length(dRef[(dRef[,'top']==as.numeric(tRef[t])),2])>=3) {
			lm <- mean(dRef[(dRef[,'top']==as.numeric(tRef[t])),2], na.rm=T)
			points(lm,tRef[t],cex=1.5, pch=5)
			ls <- sd (dRef[(dRef[,'top']==as.numeric(tRef[t])),2], na.rm=T)
			points(lm-ls,tRef[t],cex=1, pch='|')			
			points(lm+ls,tRef[t],cex=1, pch='I')
			lines (c(lm+ls,lm-ls),c(tRef[t],tRef[t]))
		} else {
			lines (c(fmin,fmax),c(tRef[t],tRef[t]))
		}
	}
}


plotCorePolygon<-function(dRef, title, y_lab, y_top, y_bot, x_lab, x_old, x_yng, proportion, division, colours) {

	# TRAP / SET MISSING VALUES
	if (missing(y_lab))
		y_lab <- 'core depth (cm)'
	
	if (missing(y_top)) 
		y_top <- min(data[,1])

	if (missing(y_bot)) 
		y_bot <- max(data[,1])

	if (missing(x_lab))
		x_lab <- ''

	if (missing(x_old))
		x_old <- 0

#	if (missing(x_yng))
#		x_yng <- 1

	if (missing(proportion))
		proportion <- TRUE

	if (missing(division))
		division <- 0

	if (missing(colours))
		colours <- c('red','blue','green')

	# SET CONSTANTS / PREFERENCES
	par(las=1, xaxs="i", yaxs="i", cex.axis=1.6, cex.main=1.6, cex.lab=1.6)	

	tRef <- levels(as.factor(dRef[,'top']))

	if (proportion) {
		if (x_lab=='x')
			x_lab <- 'proportion of layer'
		dRef[,2] <- dRef[,2]/dRef[,3]
	}

	if (missing(x_yng))
		x_yng <- max(dRef[,2], na.rm=T)*1.1

	# MEANS / STANDARD DEVIATIONS
	mgb <- mean(dRef[(dRef[,'top']>division),2],na.rm=TRUE)
	sgb <- sd(dRef[(dRef[,'top']>division),2],na.rm=TRUE)
	if (division > 0) {
		mlb <- mean(dRef[(dRef[,'top']<division),2],na.rm=TRUE)
		slb <- sd(dRef[(dRef[,'top']<division),2],na.rm=TRUE)
	}

	# PLOT AXES
	plot (dRef[,2],dRef[,'top'], ylim=c(y_bot,y_top), ylab=y_lab, xlim=c(x_old,x_yng), xlab=x_lab, main=paste(title), type='n')
	
	# PLOT TOP LAYER - STANDARD DEVIATIONS + MEAN
	polygon(c(mgb+sgb*2,mgb-sgb*2,mgb-sgb*2,mgb+sgb*2),c(division,division,maxDepth,maxDepth),lty=0, col=colours[3])
	polygon(c(mgb+sgb,mgb-sgb,mgb-sgb,mgb+sgb),c(division,division,maxDepth,maxDepth),lty=0, col=colours[2])
	lines(c(mgb,mgb),c(division,maxDepth),lwd=3, lty=1, col=colours[1])

	# PLOT BOTTOM LAYER - STANDARD DEVIATIONS + MEAN
	if (division > 0) {
		polygon(c(mlb+slb*2,mlb-slb*2,mlb-slb*2,mlb+slb*2),c(0,0,division,division),lty=0, col=colours[3])
		polygon(c(mlb+slb,mlb-slb,mlb-slb,mlb+slb),c(0,0,division,division),lty=0, col=colours[2])
		lines(c(mlb,mlb),c(0,division),lwd=3, lty=1, col=colours[1])
	}
	
	# PLOT POINTS + ERROR BARS
	points (dRef[,2],dRef[,'top'], cex=0.5)
	for (t in 1:length(tRef)) {
		fmin <- min(dRef[(dRef[,'top']==as.numeric(tRef[t])),2], na.rm=T)
		fmax <- max(dRef[(dRef[,'top']==tRef[t]),2], na.rm=T)
		if (length(dRef[(dRef[,'top']==as.numeric(tRef[t])),2])>=3) {
			lm <- mean(dRef[(dRef[,'top']==as.numeric(tRef[t])),2], na.rm=T)
			points(lm,tRef[t],cex=1.5, pch=5)
			ls <- sd (dRef[(dRef[,'top']==as.numeric(tRef[t])),2], na.rm=T)
			points(lm-ls,tRef[t],cex=1, pch=3)			
			points(lm+ls,tRef[t],cex=1, pch=3)
			lines (c(lm+ls,lm-ls),c(tRef[t],tRef[t]))
		} else {
			lines (c(fmin,fmax),c(tRef[t],tRef[t]))
		}
	}
}

# disOffset = distance off actual value in y units
# boxWidth = width of box in y units
# pns = (logical) plot number of samples in layer
######################
plotCoreLevelsError <- function(dataMatrix,colXaxis,colYaxis,disOffset,boxWidth,colour,colPoly,colOutline,pns) {

	if (missing(boxWidth))
		boxWidth <- 1

	if (missing(disOffset))
		disOffset <- 0

	if (missing(colPoly))
		colPoly <- 'grey50'

	if (missing(colOutline))
		colOutline <- 'white'

	if (missing(colour))
		colour <- 'black'

	if (missing(pns))
		pns<- FALSE

	for (i in levels(as.factor(dataMatrix[,colYaxis]))) {
	
		a <- quantile(dataMatrix[(dataMatrix[,colYaxis]==i),colXaxis],c(0.975,0.75,0.5,0.25,0.025), na.rm=T)
	
		y <- as.numeric(i) + disOffset
		
		if (pns)
			text (0,y,length(dataMatrix[(dataMatrix[,colYaxis]==i),colXaxis]))
		
		polygon(c(a[2],a[4],a[4],a[2]),c(y+boxWidth,y+boxWidth,y-boxWidth,y-boxWidth), border=colour, col=colPoly)
		lines(c(a[1],a[5]),c(y,y), col=colour)
		lines(c(a[1],a[1]),c(y+boxWidth,y-boxWidth), col=colour)
		lines(c(a[5],a[5]),c(y+boxWidth,y-boxWidth), col=colour)
		lines(c(a[3],a[3]),c(y+(1.5*boxWidth),y-(1.5*boxWidth)), lwd=3, col=colour)
	}

	points(dataMatrix[,colXaxis],dataMatrix[,colYaxis]+disOffset, pch=21, bg=colour, lwd=0.5, cex=1.0, col='white')
	points(inferedAges[,agesToUse[2]], inferedAges[,'top'],col=colOutline, bg=colLine, pch=21, cex=1.0, lwd=0.5)

}