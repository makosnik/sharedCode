#	Function for plotting cores...
#		vertical plot with area fill.
#		"jagged" layers rather than smooth points/lines.
#		- 1st column with core depth, remaining columns with data to be plotted.
###############################################################################################
plotAreaVertical<-function(data, title, y_lab, y_top, y_bot, x_lab, x_old, x_yng, colors, proportion) {

	# TRAP / SET MISSING VALUES
	if (missing(proportion))
		proportion <- TRUE

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

	if (missing(x_yng))
		if (proportion)
			x_yng <- 1
		else
			x_yng <- 100
		
	if (missing(colors))
		colors <- c('black','grey20','grey40','grey60','grey80','white','green')

	# SET CONSTANTS / PREFERENCES
	par(las=1, xaxs="i", yaxs="i", cex.axis=1.6, cex.main=1.6, cex.lab=1.6)	

	if (proportion)
		if (x_lab=='x')
			x_lab <- 'proportion of layer'

	# THE FOLLOWING FAILS IF ONLY 1 ROW IS IN THE DATASET

	# DETERMINE SIZE OF DATASET
	numCols <- length(data[1,])
	numRows <- length(data[,1])

	# MAKE NEW DATA MATRIX 
	#	with a row added at each end for polygon completion
	#	2 x as many datapoints so the area is correct for the entire layer
	data2 <- matrix(data=0,ncol=numCols,nrow=numRows*2+2)
	colnames(data2) <- colnames(data)

	# determine out which values are null, sub 0
	nar<-0
	naRows <- vector(mode='numeric', length=numRows)
	for (i in 1:numCols)
		for (j in 1:numRows)
			if (is.na(data[j,i])) {
				data[j,i] <-0
				if (i==1) {
					nar <- nar+1
					naRows[nar] <- j
				}
			}

	# change data to culumlative distribution.
	for (i in 3:numCols)
		data[,i] <- data[,i]+data[,i-1]

	# change data to proportion (if needed).
	if (proportion)
		for (i in 2:numCols)
			data[,i] <- data[,i]/data[,length(data[1,])]

	# load data into polygon data matrix
	j <- 1
	for (i in 1:(numRows)) {
		j <- j+1
		data2[j,] <- data[i,]
		j <- j+1
		data2[j,] <- data[i,]
	}
	data2[1,] <- 0
#	data2[1,1] <- min(data[1,1])
	data2[1,1] <-0
	data2[numRows*2+2,] <- 0
#	data2[numRows*2+2,1] <- max(data[,1]+(data[length(data[,1]),1]-data[(length(data[,1])-1),1]))
	data2[numRows*2+2,1] <- max(data[,1])
	
	# shift the core depths by 1 row so that the "stair steps" work as intended.
#	for (i in 1:(length(data2[,1])-1))
#		data2[i,1] <- data2[i+1,1]
	for (i in (length(data2[,1])-1):2)
		data2[i,1] <- data2[i-1,1]
	
	# PLOT DATA (Finally!)
	plot(1:1,type="n", ylim=c(y_bot,y_top), ylab=y_lab, xlim=c(x_old,x_yng), xlab=x_lab, main=title)
	
	for (i in (numCols):2) {
		polygon (data2[,i],data2[,1], col=colors[i], lty=0)
		lines (data2[,i],data2[,1], col='black', lty=1)
	}
		
}