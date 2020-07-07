##
## PLOTS TWO LINES ON PHANEROZOIC AXES.
##

timeplot<-function(dataRight, dataLeft, dataTime, title, y_lab, y_top, y_bot, x_old, x_yng, colRight, colLeft, file_out, plotPoints, plotLines) {
	
	source("/Volumes/shared/svnSifr/analyses/sharedCode/PhanerozoicPlotBox.R")
	
	if (missing(dataTime)) 
		dataTime <- binMids

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
	
	if (!missing(dataLeft)) 
	if (missing(colLeft))
		colLeft <- c('black','grey75','grey50','grey25','pink')
	if (!missing(dataRight)) 
	if (missing(colRight))
		colRight <- rainbow(length(dataRight[1,]))
	
	phanerozoicPlotBox(y_top, y_bot, x_old, x_yng, y_lab)

	if (! missing(dataRight))
		axis(4,at=c(0,y_scale*.25,y_scale*.5,y_scale*.75,y_scale),labels=c('0%','25%','50%','75%','100%'), col=colRight[1], col.lab=colRight[1], col.axis=colRight[1], cex.axis=1.6, cex.lab=1.6)
		
	
	if (! missing(dataRight))
		for (i in 1:length(dataRight[1,])) {
			if (plotLines)
				lines(dataTime,dataRight[,i]*y_scale, lty=1, col=colRight[i])
			if (plotPoints)
				points(dataTime,dataRight[,i]*y_scale, pch=18, col=colRight[i])
		}
		
	if (! missing(dataLeft))
		for (i in 1:length(dataLeft[1,])) {
			if (plotLines)
				lines(dataTime,dataLeft[,i], lty=1, col=colLeft[i])
			if (plotPoints)
				points(dataTime,dataLeft[,i], pch=19, col=colLeft[i])
		}
	
	
#	legend()
	
	if (file_out == TRUE) {
		dev.off()
	}
}