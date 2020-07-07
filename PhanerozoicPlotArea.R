timeplotarea<-function(dataArea, dataTime, title, y_lab, y_top, y_bot, x_old, x_yng, colors) {

	source("/Volumes/shared/svnSifr/analyses/sharedCode/PhanerozoicPlotBox.R")

	# MESS WITH DATA TO GET POLYGONS
	numCols <- length(dataArea[1,])
	numRows <- length(dataArea[,1])

	# new data matrix with a row added at each end for polygon completion
	binMid2 <- vector(mode='numeric',length=numRows+2)
	data2 <- matrix(data=0,ncol=numCols,nrow=numRows+2)
	colnames(data2) <- colnames(dataArea)

	# figure out which time bins are null, sub 0
	nar<-0
	naRows <- vector(mode='numeric', length=numRows)
	for (i in 1:numCols)
		for (j in 1:numRows)
			if (is.na(dataArea[j,i])) {
				dataArea[j,i] <-0
				if (i==1) {
					nar <- nar+1
					naRows[nar] <- j
				}
			}

	#  - change data to culumlative distribution.
	for (i in 2:numCols)
		dataArea[,i] <- dataArea[,i]+dataArea[,i-1]

	# load data into polygon data matrix
	for (i in 1:(numRows)) {
		data2[i+1,] <- dataArea[i,]
		binMid2[i+1] <- dataTime[i]
	}
	data2[1,] <- 0
	data2[numRows+2,] <- 0
	binMid2[1] <- max(dataTime)
	binMid2[numRows+2] <- min(dataTime)

	if (missing(x_old)) x_old <- max(binMid2)
	if (missing(x_yng)) x_yng <- min(binMid2)

	
	if (missing(y_top)) y_top <- 0
	if (missing(y_bot)) y_bot <- 0

	if (min(data2[,1])<y_bot) y_bot <- min(data2[,1])
	if (max(data2[,numCols])>y_top) y_top <- max(data2[,numCols])


	# DEFAULTS
	blockColor <- 'white'

	if (missing(y_lab))
		y_lab <- ''

	if (missing(colors))
		colors <- c('black','grey75','grey50','grey25','pink')
		
	phanerozoicPlotBox(y_top, y_bot, x_old, x_yng)	
	axis(side=2,at=seq(y_top, y_bot, length.out=6))
	
	# - plot data
	for (i in numCols:1)
		polygon(binMid2,data2[,i], col=colors[i])

	# - block out areas w/o data
#	for (i in 1:(numRows))
#		if (naRows[i]>0)
#			polygon (c(binMid2[(naRows[i])],binMid2[(naRows[i])],binMid2[(naRows[i]+2)],binMid2[(naRows[i]+2)]),c(y_bot,y_top,y_top,y_bot), col=blockColor, lty=0)
#	polygon (c(min(binMid2),min(binMid2),x_yng,x_yng),c(y_bot,y_top,y_top,y_bot), col=blockColor, lty=0)
#	polygon (c(max(binMid2),max(binMid2),x_old,x_old),c(y_bot,y_top,y_top,y_bot), col=blockColor, lty=0)
	
#	legend(510,(((y_top-y_bot)/2)+y_bot),colnames(data2),colors, horiz=FALSE, bty='n', cex=1)
			
}