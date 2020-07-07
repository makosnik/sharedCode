#	$Revision: 607 $

IQ0 <- 0.5
SD1 <- 0.6827
SD2 <- 0.9545
SD3 <- 0.9973


qIQ0 <- c(0.50-(IQ0/2),0.5,0.5+(IQ0/2))
qSD1 <- c(0.50-(SD1/2),0.5,0.5+(SD1/2))
qSD2 <- c(0.50-(SD2/2),0.5,0.5+(SD2/2))
qSD3 <- c(0.50-(SD3/2),0.5,0.5+(SD3/2))


#######################################################################################################
## REQUIRED DATA
#select sample_no, layer_no, top, depth, fraction_no, sieve_size, sediment_dry_weight_g from fractions join layers using (layer_no) where sample_no=592 into outfile "/tmp/Careel-coreSed.txt";

# NEEDS:
# 
# DATA.FRAC - DATABASE FRACTION DUMP (SEE EXAMPLE SQL LINE ABOVE)
# PLOT.NAME - NAME FOR PLOT FILE
# PLOT.MAIN - MAIN FOR PLOTS

coreFractions <- unique(as.character(DATA.FRAC[,'sieve_size']))
coreLayers <- unique(DATA.FRAC[,'top'])

columnNames <- c('top','depth',coreFractions,'total')
coreSed <- matrix(nrow=length(coreLayers),ncol=length(columnNames))
coreSed[,1] <- coreLayers
colnames(coreSed) <- columnNames

for (i in 1:length(coreLayers)) {
	coreSed[i,'depth'] <- DATA.FRAC[((DATA.FRAC[,'sieve_size']=='4') & (DATA.FRAC[,'top']== coreLayers[i])),'depth']	
	for (s in coreFractions) {
		coreSed[i,as.character(s)] <- DATA.FRAC[((DATA.FRAC[,'sieve_size']==s) & (DATA.FRAC[,'top']== coreLayers[i])),'sediment_dry_weight_g']
	}
}

coreSed[,'total'] <- rowSums(coreSed[,as.character(coreFractions)],na.rm=TRUE,dims=1)
coreSed[(coreSed[,'total'] == 0),'total'] <- NA

coreSedPro <- coreSed
for (s in coreFractions)
	coreSedPro[,s] <- coreSed[,s]/coreSed[,'total']
proMin <- rep(0,length= length(coreLayers))
proMax <- rep(1,length= length(coreLayers))
coreSedPro <- cbind(coreSedPro,proMax,proMin)
coreFractionsPro <- c('proMax', coreFractions, 'proMin')
sedCumPro <- coreSedPro
for (f in seq(length(coreFractions)-1,1,by=-1)) {
	sedCumPro[,coreFractions[f]] <- sedCumPro[,coreFractions[f]]+ sedCumPro[,coreFractions[f+1]]
}

#######################################################################################################
pdf(paste(PLOT.NAME,sep='/'), width=pageWidthTwo, height=pageHeight/2)
colLine <- 'black'
colText <- 'black'
cexDefault <- 1.0
par (cex=cexDefault, fg=colLine, mar=c(4,4,1,1), mgp=c(2.2,0.7,0))

par(mfcol=c(1,2))

ytop <- coreSed[1,'top']
ybot <- coreSed[nrow(coreSed),'top']+coreSed[nrow(coreSed),'depth']+10

ytop <- min(coreSed[,'top'])-2
ybot <- max(coreSed[,'top'])+20

xMin <- min(coreSed[,coreFractions]/coreSed[,'depth'], na.rm=TRUE)*0.95
xMax <- max(coreSed[,'total']/coreSed[,'depth'], na.rm=TRUE)*1.05

if (xMin < 0.025)
	xMin <- 0.025

##
##	MASS PER CM DEPTH
#######################################################################################################

plot(coreSed[,'total']/coreSed[,'depth'],coreSed[,'top']+coreSed[,'depth']/2, type='l', main=PLOT.MAIN, log='x', xlim=c(xMin,xMax), xlab='Sediment mass (g/cm)', ylab='Excavation depth (m)', yaxt='n', yaxs='i', ylim=c(ybot,ytop), lwd=2, pch=20)
axis(2,at=seq(0,160, by=20),labels=formatC(seq(0,1.6, by=0.20),format='f',digits=1),col=colLine, col.axis=colText, las=1)
lineType <- 2
points(coreSed[,'total']/coreSed[,'depth'],coreSed[,'top']+coreSed[,'depth']/2,cex=1.0, pch=20)

pointType <- 21
for (s in coreFractions) {
	points(coreSed[,s]/coreSed[,'depth'],coreSed[,'top']+coreSed[,'depth']/2,cex=0.5, pch=pointType, bg='grey')
	lines(coreSed[,s]/coreSed[,'depth'],coreSed[,'top']+coreSed[,'depth']/2,lty=lineType, lwd=0.5)
	lineType <- lineType + 1
	pointType <- pointType + 1
}

legend(xMin,ybot*0.6,legend=c('Total',paste('>',coreFractions,sep='')), lty=seq(1,length(coreFractions)+1, by=1), lwd=c(2,rep(0.5,length(coreFractions))), pch=seq(from=20, length.out=length(coreFractions)+1, by=1), pt.cex=c(1,rep(0.5,length(coreFractions))), pt.bg='grey')

oneProMaxTop <- coreSedPro[1,'1']+coreSedPro[1,'2']/2
oneProMaxBot <- coreSedPro[nrow(coreSedPro),'1']+coreSedPro[nrow(coreSedPro),'2']/2

##
##	PROPORTION MASS PER LAYER
#######################################################################################################

plot(sedCumPro[,'1'], sedCumPro[,'top']+ sedCumPro[,'depth'],ylim=c(ybot,ytop), xlim=c(0,1), type='n', xaxs='i', yaxs='i', ylab='',xlab='Proportion dry mass', main=PLOT.MAIN, col=colLine, col.axis=colLine,col.sub=colLine,col.main=colText,col.lab=colText, yaxt='n')
axis(2,at=seq(0,160, by=20),labels=formatC(seq(0,1.6, by=0.20),format='f',digits=1),col=colLine, col.axis=colText, las=1)

cfp <- coreFractionsPro
colours <- gray(0:length(coreFractionsPro) / length(coreFractionsPro))

#c('grey20','grey30','grey50','grey70','grey80')

## fill areas
for (i in 1:nrow(sedCumPro)) {
	for (f in 2:length(cfp)-1) {
		colFill <- colours[f+1]
		colBord <- colours[f-1]
		polygon(c(sedCumPro[i, cfp[f]], sedCumPro[i, cfp[f+1]], sedCumPro[i, cfp[f+1]], sedCumPro[i, cfp [f]]),c(sedCumPro[i,'top'], sedCumPro[i,'top'], sedCumPro[i,'top']+ sedCumPro[i,'depth'], sedCumPro[i,'top']+ sedCumPro[i,'depth']),col=colFill, border=colBord)
	}
}

# layer lines
for (i in 1:nrow(coreSedPro))
	lines(c(0,1),c(coreSedPro[i,'top'], coreSedPro[i,'top']), lty=3)
lines(c(0,1),c(coreSedPro[nrow(coreSedPro),'top']+coreSedPro[nrow(coreSedPro),'depth'], coreSedPro[nrow(coreSedPro),'top']+coreSedPro[nrow(coreSedPro),'depth']), lty=3)


## ADD MEDIAN PROPORTION (VERTICAL LINES)
previousLineValue <- 0
for (f in as.character(sort(as.numeric(coreFractions)))) {
	lineValue <- median(sedCumPro[,f], na.rm=TRUE)
	proValue <- median(coreSedPro[,f], na.rm=TRUE)
	lines(c(lineValue, lineValue),c(coreSed[nrow(coreSed),'top']+coreSed[nrow(coreSed),'depth']+2.5,coreSed[1,'top']), lty=2, lwd=2, col=colLine)
	text(((lineValue-previousLineValue)/2)+ previousLineValue,coreSed[nrow(coreSed),'top']+coreSed[nrow(coreSed),'depth']+12,paste(">",f,sep=''),pos=3, font=4, col=colText, cex=0.7)
	text(((lineValue-previousLineValue)/2)+ previousLineValue,coreSed[nrow(coreSed),'top']+coreSed[nrow(coreSed),'depth']+18,round(proValue,2),pos=3, font=4, col=colText, cex=0.7)
	previousLineValue <- lineValue
}

dev.off()