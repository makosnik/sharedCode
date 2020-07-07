
plotPolygon <- function(plotData,timeData,fillColor,lineWeight) {
		
	yData <- vector(mode='numeric',length=(2*length(plotData[,1])))
	xData <- vector(mode='numeric',length=(2*length(timeData)))
	c <- 0		

	for (f in 1:length(plotData[,1])) {
#		if (!is.na(plotData[f,1])) {
			c <- c+1
			yData[c] <- plotData[f,1]
			xData[c] <- timeData[f]
#		}
	}

	for (f in length(plotData[,2]):1) {
#		if (!is.na(plotData[f,2])) {
			c <- c+1
			yData[c] <- plotData[f,2]
			xData[c] <- timeData[f]
#		}
	}

	xData <- xData[!is.na(yData)]
	yData <- yData[!is.na(yData)]

	
#	yData <- subset(yData,xData>=0)		
#	xData <- subset(xData,xData>=0)
	
	polygon(xData,yData,col=fillColor,lwd=lineWeight)

}

#
# PLOTS MATRIX OF LINES ON PHANEROZOIC AXIS
#
timeplot3<-function(plotData, dataTime, title, y_lab, y_top, y_bot, x_old, x_yng, Data_col, file_out, plotPoints, plotLines) {
		
	source("/Volumes/shared/svnSifr/analyses/sharedCode/PhanerozoicPlotBox.R")

	if (missing(title)) 
		title <- ''

	if (missing(plotPoints)) 
		plotPoints <- FALSE

	if (missing(plotLines)) 
		plotLines <- TRUE
	
	if (missing(file_out)) 
		file_out <- TRUE
	if (file_out == TRUE) {
		ps.options(paper="letter", horizontal=FALSE, onefile=FALSE, family="Helvetica", pointsize=10, encoding="ISOLatin1", width=8, height=8)
	
		postscript(paste('tp-',title,".eps",sep=''))
		par(las=1, pin=c(4,3))
	}

	if (missing(Data_col))
		Data_col <- 'black'
	

	# Color forground
	colF <- 'black'
	# Color background
	colB <- 'white'

	colCI <- 'grey50'
	colSE <- 'grey25'

	
	phanerozoicPlotBox(y_top, y_bot, x_old, x_yng)
	
	if (plotLines) {
		
		if (length(plotData[1,]) >= 3)
			plotPolygon(plotData[,c(1,(length(plotData[1,])-0))],dataTime,'grey75',0.1)
		if (length(plotData[1,]) >= 5)
			plotPolygon(plotData[,c(2,(length(plotData[1,])-1))],dataTime,'grey50',0.1)
		if (length(plotData[1,]) >= 7)
			plotPolygon(plotData[,c(3,(length(plotData[1,])-2))],dataTime,'grey25',0.1)

		lines(dataTime,plotData[,ceiling(length(plotData[1,])/2)], lty=1, lwd=1.5, col=Data_col)
	}
	
	if (plotPoints) {
#		points(dataTime,plotData[,2], pch=19, cex=0.75, col=colSE)
#		points(dataTime,plotData[,4], pch=19, cex=0.75, col=colSE)

		points(dataTime,plotData[,1], pch=25, cex=0.75, col=colCI, bg=colCI)
		points(dataTime,plotData[,3], pch=24, cex=0.75, col=colCI, bg=colCI)

		points(dataTime,plotData[,2], pch=19, cex=1.00, col=colF)
	
	}
		
	if (file_out == TRUE) {
		dev.off()
	}
}

#
# PLOTS MATRIX OF LINES ON PHANEROZOIC AXIS
#
timeplot2<-function(plotData, dataTime, title, y_lab, y_top, y_bot, x_old, x_yng, Data_col, file_out, plotPoints, plotLines) {
		
	source("/Volumes/shared/svnSifr/analyses/sharedCode/PhanerozoicPlotBox.R")

	if (missing(title)) 
		title <- ''

	if (missing(plotPoints)) 
		plotPoints <- FALSE

	if (missing(plotLines)) 
		plotLines <- TRUE
	
	if (missing(file_out)) 
		file_out <- TRUE
	if (file_out == TRUE) {
		ps.options(paper="letter", horizontal=FALSE, onefile=FALSE, family="Helvetica", pointsize=10, encoding="ISOLatin1", width=8, height=8)
	
		postscript(paste('tp-',title,".eps",sep=''))
		par(las=1, pin=c(4,3))
	}

	if (missing(Data_col))
		Data_col <- 'black'
	
	# Color forground
	colF <- 'black'
	# Color background
	colB <- 'white'

	colCI <- 'grey50'
	colSE <- 'grey25'

	phanerozoicPlotBox(y_top, y_bot, x_old, x_yng)
	
	if (plotLines) {
#		lines(dataTime,plotData[,2], lty=1, lwd=1, col=colSE)
#		lines(dataTime,plotData[,4], lty=1, lwd=1, col=colSE)

		lines(dataTime,plotData[,1], lty=1, lwd=1, col=colCI)
		lines(dataTime,plotData[,3], lty=1, lwd=1, col=colCI)

		lines(dataTime,plotData[,2], lty=1, lwd=1.5, col=Data_col)
	}
	
	if (plotPoints) {
#		points(dataTime,plotData[,2], pch=19, cex=0.75, col=colSE)
#		points(dataTime,plotData[,4], pch=19, cex=0.75, col=colSE)

		points(dataTime,plotData[,1], pch=25, cex=0.75, col=colCI, bg=colCI)
		points(dataTime,plotData[,3], pch=24, cex=0.75, col=colCI, bg=colCI)

		points(dataTime,plotData[,2], pch=19, cex=1.00, col=colF)
	
	}
		
	if (file_out == TRUE) {
		dev.off()
	}
}
