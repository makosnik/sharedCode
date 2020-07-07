## DOWN-CORE TRENDS IN EACH TAXON - AGE
###################################################
source('../sharedCode/plotCore.R')
boxWidth <- 1
plotSize <- 5
pSpacer <- 0.98		# relative position on the x-axis of the p-value plots (1=left, 0=right)

#xColumn <- 'wmAll'
xLabel <- 'mean age (years)'
#yColumn <- 'plotDepth'

xmax <- 1
xmin <- 0
xmax <- (ceiling(max(dataToUse[, xColumn], na.rm=TRUE)/100)*100)+20
xmin <- (floor(min(dataToUse[, xColumn], na.rm=TRUE)/100)*100)-20
ymax <- (ceiling(max(dataToUse[, yColumn], na.rm=TRUE)/10)*10)+5
ymin <- -5

# SET PAGE HORIZONTAL
# ps.options(paper="A4", horizontal=TRUE, onefile=TRUE, family="Helvetica", pointsize=14, encoding="ISOLatin1.enc")
 ps.options(paper="A4", horizontal=FALSE, onefile=TRUE,family= fontFamily, pointsize=12, encoding="ISOLatin1.enc", width=pageWidthTwo, height=pageHeight/2)

for (core in SITES) {
	
	# DATA BY CORE
	dataCore <- subset(dataToUse, sample_no==core)

	# OPEN PLOT FILE
	postscript(paste(FILEOUT,"coreAge-",core,".eps",sep=''))
	par(mfrow=c(1,4), mar=c(3,3,2,0), las=0, mgp=c(2,1,0))
	
	# PLOT FORMAT
	rows <- 1
#	cols <- length(TAXA)
	cols <- 2
	nplots <- rows * cols
#	layout(matrix(seq(nplots), rows, cols, byrow=TRUE), width=rep(lcm(plotSize),cols), height=rep(lcm(plotSize*2),rows), respect=TRUE)
	

	# INITIALIZE FIGURE LETTERS
	let <- 0

#	TAXA <- levels(factor(dataToUse[,'material', drop=TRUE]))
	for (taxon in TAXA) {

		# DATA BY TAXON
		dataTaxon <- subset(dataCore, taxon==taxon)
		
		# IF DATA...
		if (nrow(dataTaxon) > 1) {
		
			# MAKE PLOT TITLE
			let <- let+1
			
			if (let == 1) {
				yLabel <- 'core depth (cm)'
			} else {
				yLabel <- ''
			}
			
			# MAKE PLOT
			plot(dataTaxon[,xColumn],dataTaxon[,yColumn], ylim=c(ymax,ymin), xlim=c(xmax,xmin), ylab=yLabel, xlab=xLabel, type='n', xaxs='i', yaxs='i', xaxt='n', main='', pch=20)

			# TITLE
			title(paste(LETTERS[let],'.',sep=''),adj=0)
#			title(main=taxon, font.main=4)
			title (main=taxon, adj=0.5, font.main=3)
			
			# X-AXIS
#			axis(1,c(0,'',2000,'',4000,''), at=seq(0,5000,1000))
			axis(1,seq(0,1200,200), at=seq(0,1200,200))
			
			# ADD ERROR BARS
			plotCoreLevelsError(dataTaxon, xColumn, yColumn, pns=FALSE, boxWidth=2)

			# ADD DIVDER LINE @ 25CM
			lines(c(xmax,xmin),c(25,25), lty=2, lwd=1)

			# CALC/PRINT - KRUSKAL TEST - ALL CORE LAYERS
			if (min(dataTaxon[, yColumn])<25) {
				pk <- kruskal.test(dataTaxon[, xColumn], dataTaxon[,yColumn])
				if (pk$p.value >= 0.01) {
					text(xmax*pSpacer,5, paste('p =',round(pk$p.value,2)), pos=4)
				} else {
					text(xmax*pSpacer,5, 'p < 0.01', pos=4)
				}
			}

			# CALC/PRINT - KRUSKAL TEST - LAYERS DEEPER THAN 25CM
			pk <- kruskal.test(dataTaxon[(dataTaxon[,yColumn]>25), xColumn], dataTaxon[(dataTaxon[,yColumn]>25), yColumn])
			if (pk$p.value >= 0.01) {
				text(xmax*pSpacer,30, paste('p =',round(pk$p.value,2)), pos=4)
			} else {
				text(xmax*pSpacer,30, 'p < 0.01', pos=4)			
			}
		}
	}
	
	# CLOSE PLOT FILE
	dev.off()
}
