################################  FILE LICENSE  ################################
#
#	SVN repository at sifr.science.mq.edu.au
#	$Revision: 1040 $
#
#	This file is copyright (C) 2014 Matthew Kosnik
#
#	This program is free software; you can redistribute it and/or modify it 
#	under the terms of version 2 the GNU General Public License as published 
#	by the Free Software Foundation.
#
#	This program is distributed in the hope that it will be useful, but WITHOUT
#	ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
#	FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for 
#	more details.
#
#	To view a copy of the license go to:
#	http://www.fsf.org/copyleft/gpl.html
#	To receive a copy of the GNU General Public License write the Free Software
# 	Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
################################################################################

## NEEDS TO BE SET:

# PATH

# FILE.DAT

# taxon






##	INPUT FILE - FROM DATABASE
# SELECT layers.*, fractions.*, specimens.*,measure2.*,c14Dates.*,aarPeakAreas.* from specimens left join measure2 using (specimen_no) left join c14Dates using (specimen_no) left join aarPeakAreas using (specimen_no) left join fractions using (fraction_no) left join layers using (layer_no) 
# WHERE taxon_no=17469 into outfile '/tmp/sydPeronella.txt';
# WHERE taxon_no=10326 into outfile '/tmp/sydFulvia.txt';
# WHERE taxon_no=17464 into outfile '/tmp/sydTimoclea.txt';
#layer_no	sample_no	top	depth	bag_num	collector_no	date_collected	mass_gm	collectors	comments	fraction_no	layer_no	sieve_size	total_dry_weight_g	bag_weight_g	sediment_dry_weight_g	tmp_container	date_picked	comments	use_it	full	specimen_no	user_no	taxon_no	taxon_reso	fraction_no	number	skip_in_count	specimen_type	vital	completeness	breakage	breakage_time	condition	quality	age	encrustation	comments	specimen_fate	dom	measure2_no	user_no	specimen_no	date_measured	measured_object	x_mm	y_mm	z_mm	thick_mm	sinus_x_mm	sinus_y_mm	mass_mg	density	divide_mass_by	bias	dom	c14date_no	specimen_no	labCode	ansto_num	d13c	pMC	e1Sigma	d14C	ed14C	cAgeBP	cAgeSigma	2sYng	2sOld	2sMedian	ageMedian	ageYng	ageOld	resAge	calCurve	date_prep	comments	area_no	specimen_no	ual_num	material	location	treatment	weight	flag	gain	L_Asp	D_Asp	L_Glu	D_Glu	L_Ser	D_Ser	Gly	L_Ala	hArgL	D_Ala	L_Val	D_Val	L_Phe	L_Ile	D_Phe	L_Leu	D_Alle	D_Leu	hArgL_source


# SET PATH TO LOCATION OF FILES
setwd('/Volumes/shared/svnSifr/analyses/trunk/papers/')
PATH <- './2014sydDates'

# THE OUTPUT FROM THE ANALYSES WILL BE PUT IN THIS DIRECTORY
# can be specified directly or relative to the project directory as below 
PATH.OUT <- paste(PATH,'output',sep='/')
if(!paste(PATH.OUT) %in% list.dirs(PATH)) dir.create(PATH.OUT)

# PATH CONTAINING DATA FILES
PATH.DAT <- paste(PATH,'data',sep='/')

# FILE CONTAINING DATA
# can be specified directly or relative to the project directory as below
# if you don't have a *.csv file with a single header row then you will need tweak the loading line too
FILE.DAT <- paste(PATH.DAT,'sydFulvia.txt',sep='/')

## LOAD DATA
#DATA.RAW <- read.csv(FILE.DAT, header=TRUE, na.strings=c("\\N","NA"))
DATA.RAW <- read.table(FILE.DAT, header=TRUE, na.strings=c("\\N","NA"), sep='\t', quote='')
areaRaw <- DATA.RAW

##	SUBSET FILE IF NEEDED:

## ONLY USE top edge samples
taxon <- 'Peronella peronii'
taxon <- 'Timoclea (Chioneryx) cardioides'
taxon <- 'Fulvia tenuicostata'
#areaRaw <- areaRaw[(areaRaw[,'taxon_no']=='17469'),]
areaRaw <- areaRaw[!is.na(areaRaw[,'D_Asp']),]

#areaRaw <- areaRaw[is.na(areaRaw[,'location']) | (areaRaw[,'location']=='top edge'),]


TAXA <- unique(as.factor(areaRaw[,'taxon_no']))
TOPS <- unique(as.factor(areaRaw[,'top']))
LAYERS <- unique(as.factor(areaRaw[,'layer_no']))
FRACTIONS <- unique(as.factor(areaRaw[,'fraction_no']))



###################################################
# SET - KEY PARAMETERS / CONSTANTS
###################################################
iniA <- c('Asp','Glu','Ser','Ala','Val','Phe','Leu','A/I')
namA <- c('Aspartic','Glutamic','Serine','Alanine','Valine','Pheylalanine','Leucine','Alle/Ile')
iniA2 <- c('Asp','Glu','Ser','Ala','Val','Phe','Leu','AI')

colD <- paste('cD', iniA2,sep='_')
colL <- paste('cL', iniA2,sep='_')
colR <- paste(iniA2,'DL', sep='_')

SD1 <- 0.6827
SD2 <- 0.9545
qSD1 <- c(0.50-(SD1/2),0.5,0.5+(SD1/2))
qSD2 <- c(0.50-(SD2/2),0.5,0.5+(SD2/2))
qSD <- c(0.50-(SD2/2),0.50-(SD1/2),0.5,0.5+(SD1/2),0.5+(SD2/2))

###################################################
#	CHECK SIZES
###################################################

areaRaw[((!is.na(areaRaw[,'z_mm'])) & (!is.na(areaRaw[,'measured_object'])) & (areaRaw[,'measured_object']=='whole')),'z_mm'] <- areaRaw[((!is.na(areaRaw[,'z_mm'])) & (!is.na(areaRaw[,'measured_object'])) & (areaRaw[,'measured_object']=='whole')),'z_mm'] / 2


gm_size <- (areaRaw[,'x_mm']*areaRaw[,'y_mm']*areaRaw[,'z_mm'])^(1/3)

areaRaw <- cbind(areaRaw,gm_size)

areaRaw[,'mass_mg'] <-  areaRaw[,'mass_mg']/areaRaw[,'divide_mass_by']
areaRaw[,'divide_mass_by'] <-  1

par(mfcol=c(2,2))

plot(areaRaw[,'x_mm'],areaRaw[,'y_mm'])
plot(areaRaw[,'x_mm'],areaRaw[,'z_mm'])
plot(areaRaw[,'y_mm'],areaRaw[,'z_mm'])
plot(areaRaw[,'mass_mg'],areaRaw[,'gm_size'], log='xy')

a <- lm(log(areaRaw[,'mass_mg']/areaRaw[,'divide_mass_by']) ~ log(areaRaw[,'gm_size']))
plot(a)

###################################################
#	CALCULATE - D & L CONCENTRATIONS
###################################################
cD_Asp <- (areaRaw[,'D_Asp']/areaRaw[,'hArgL'])*200
cL_Asp <- (areaRaw[,'L_Asp']/areaRaw[,'hArgL'])*200
cD_Glu <- (areaRaw[,'D_Glu']/areaRaw[,'hArgL'])*200
cL_Glu <- (areaRaw[,'L_Glu']/areaRaw[,'hArgL'])*200
cD_Ser <- (areaRaw[,'D_Ser']/areaRaw[,'hArgL'])*200
cL_Ser <- (areaRaw[,'L_Ser']/areaRaw[,'hArgL'])*200
cD_Ala <- (areaRaw[,'D_Ala']/areaRaw[,'hArgL'])*200
cL_Ala <- (areaRaw[,'L_Ala']/areaRaw[,'hArgL'])*200
cD_Val <- (areaRaw[,'D_Val']/areaRaw[,'hArgL'])*200
cL_Val <- (areaRaw[,'L_Val']/areaRaw[,'hArgL'])*200
cD_Phe <- (areaRaw[,'D_Phe']/areaRaw[,'hArgL'])*200
cL_Phe <- (areaRaw[,'L_Phe']/areaRaw[,'hArgL'])*200
cD_Leu <- (areaRaw[,'D_Leu']/areaRaw[,'hArgL'])*200
cL_Leu <- (areaRaw[,'L_Leu']/areaRaw[,'hArgL'])*200
cD_AI  <- (areaRaw[,'D_Alle']/areaRaw[,'hArgL'])*200
cL_AI  <- (areaRaw[,'L_Ile']/areaRaw[,'hArgL'])*200

###################################################
#	CALCULATE - D/L RATIOS
###################################################
Asp_DL <- areaRaw[,'D_Asp']/areaRaw[,'L_Asp']
Glu_DL <- areaRaw[,'D_Glu']/areaRaw[,'L_Glu']
Ser_DL <- areaRaw[,'D_Ser']/areaRaw[,'L_Ser']
Ala_DL <- areaRaw[,'D_Ala']/areaRaw[,'L_Ala']
Val_DL <- areaRaw[,'D_Val']/areaRaw[,'L_Val']
Phe_DL <- areaRaw[,'D_Phe']/areaRaw[,'L_Phe']
Leu_DL <- areaRaw[,'D_Leu']/areaRaw[,'L_Leu']
AI_DL <- areaRaw[,'D_Alle']/areaRaw[,'L_Ile']

areaRaw <- cbind(areaRaw,cD_Asp, cL_Asp, cD_Glu, cL_Glu, cD_Ser, cL_Ser, cD_Ala, cL_Ala, cD_Val, cL_Val, cD_Phe, cL_Phe, cD_Leu, cL_Leu, cD_AI, cL_AI)
areaRaw <- cbind(areaRaw, Asp_DL, Glu_DL, Ser_DL, Ala_DL, Val_DL, Phe_DL, Leu_DL, AI_DL)

###################################################
# DATA - SAVE CLEANED DATA
###################################################
save(areaRaw, colR, colL, colD, iniA, namA, TAXA, TOPS, file=paste(PATH.DAT,"aarAreaRaw.Rdata",sep='/'))

###

dataToUse <- areaRaw

datedSpecs <- dataToUse[!is.na(dataToUse[,'X2sMedian']),]
agesToUse <- c('X2sYng','X2sMedian','X2sOld')

#########################################################################################
##	AMINO ACID D/L CORRELATIONS
#########################################################################################
iniA3 <- c('Asp','Glu','Ala','Val','Phe','Leu','AI')
iniA3 <- c('Asp','Glu','Ala','Val','Phe')
iniA3 <- c('Asp','Glu','Ala')

colours <- c('pink','red','blue','grey')

par(mfcol=c(length(iniA3),length(iniA3)), cex=0.5, mar=c(1,1,1,1))

resid <- matrix(nrow=nrow(dataToUse),ncol=length(iniA3)*length(iniA3))
colname <- list()
i <- 0
for (x in iniA3) {
	for (y in iniA3) {
		i<-i+1
		colname[i] <- paste(x,y,sep='')
	}
}
colnames(resid) <- colname

maybe <- c(18142,18157,18177,18268,18272,18138,18131)
YesMe <- c(18137,18145,18155,18175,18161,18264,18147,18183, 18179)
NotMe <- c(18254,18150,18143)

for (x in iniA3) {
	for (y in iniA3) {
		
		if (x==y) {
			plot(1:1, type='n', xaxt='n', yaxt='n', ylab='', xlab='')
			text(1,1.1,paste(x,'D/L'), cex=3)
			text(1,0.9,taxon, cex=1.5, font=3)
			text(1,0.75,'Sydney Harbour', cex=1)
		} else {
			xLim <- c(min(dataToUse[,paste(x,'DL',sep='_')]),max(dataToUse[,paste(x,'DL',sep='_')]*1.1))
			yLim <- c(min(dataToUse[,paste(y,'DL',sep='_')]),max(dataToUse[,paste(y,'DL',sep='_')]*1.1))
			plot(dataToUse[,paste(x,'DL',sep='_')], dataToUse[,paste(y,'DL',sep='_')], xlim=xLim, ylim=yLim, xlab='',ylab='', pch=21,col='white',bg='black', lwd=0.2)
		
			a <- cor.test(dataToUse[,paste(x,'DL',sep='_')], dataToUse[,paste(y,'DL',sep='_')], method='spearman')
			b <- lm(dataToUse[,paste(y,'DL',sep='_')] ~ dataToUse[,paste(x,'DL',sep='_')])
			resid[, paste(x,y,sep='')] <- b$residuals
			abline(b)
			
			if (nrow(datedSpecs) > 0) {
				points(datedSpecs[,paste(x,'DL',sep='_')], datedSpecs[,paste(y,'DL',sep='_')],pch=19,col='blue')
				text(datedSpecs[,paste(x,'DL',sep='_')], datedSpecs[,paste(y,'DL',sep='_')], datedSpecs[,'specimen_no'],pos=4,col='blue')
			}
			
			if (length(NotMe) > 1) {
				for (c in NotMe) {
					points(dataToUse[(dataToUse[,'specimen_no']==c),paste(x,'DL',sep='_')], dataToUse[(dataToUse[,'specimen_no']==c),paste(y,'DL',sep='_')],pch=19,col='grey')
					text(dataToUse[(dataToUse[,'specimen_no']==c),paste(x,'DL',sep='_')], dataToUse[(dataToUse[,'specimen_no']==c),paste(y,'DL',sep='_')],c,pos=4,col='grey')
				}
			}

			if (length(maybe) > 1) {
				for (c in maybe) {
					points(dataToUse[(dataToUse[,'specimen_no']==c),paste(x,'DL',sep='_')], dataToUse[(dataToUse[,'specimen_no']==c),paste(y,'DL',sep='_')],pch=19,col='red')
					text(dataToUse[(dataToUse[,'specimen_no']==c),paste(x,'DL',sep='_')], dataToUse[(dataToUse[,'specimen_no']==c),paste(y,'DL',sep='_')],c,pos=4,col='red')
				}
			}

			if (length(YesMe) > 1) {
				for (c in YesMe) {
					points(dataToUse[(dataToUse[,'specimen_no']==c),paste(x,'DL',sep='_')], dataToUse[(dataToUse[,'specimen_no']==c),paste(y,'DL',sep='_')],pch=24,col='black')
					text(dataToUse[(dataToUse[,'specimen_no']==c),paste(x,'DL',sep='_')], dataToUse[(dataToUse[,'specimen_no']==c),paste(y,'DL',sep='_')],c,pos=4,col='black')
				}
			}
			
			text(max(dataToUse[,paste(x,'DL',sep='_')], na.rm=TRUE)*1,min(dataToUse[,paste(y,'DL',sep='_')], na.rm=TRUE)*1.2,paste('rho=',round(a$estimate,2)), pos=2)
			text(min(dataToUse[,paste(x,'DL',sep='_')], na.rm=TRUE)*1.2,max(dataToUse[,paste(y,'DL',sep='_')], na.rm=TRUE)*0.9,paste('r2=',round(summary(b)$r.squared,2)), pos=4)
		}
	}
}

dataToUse <- cbind(dataToUse,resid)


#########################################################################################
##	AVERAGE SPECIMEN DATA
#########################################################################################

dataOne <- dataToUse
dataOne <- dataOne[!is.na(dataOne[,'Asp_DL']),]

iniA3 <- c('Asp','Glu','Ala','Val','Phe','Leu','AI')
iniA3 <- c('Asp','Glu','Ala')
iniA3 <- c('Asp','Glu','Ala','Val','Phe')


specs <- unique(dataOne[,'specimen_no'])

for (s in specs) {
	specRows <- dataOne[(dataOne[,'specimen_no']==s),]
	if (nrow(specRows)>1) {
		for (a in iniA3) {
			specRows[1,paste(a,'DL',sep='_')] <- mean(specRows[,paste(a,'DL',sep='_')],na.rm=TRUE)
		}	
	}
	if (s == specs[1])
		dataSpec <- dataOne[1,]
	else 
		dataSpec <- rbind(dataSpec,specRows[1,])
}

dataOne <- dataSpec

#########################################################################################

maxes <- dataOne[(dataOne[,'Asp_DL']==max(dataOne[,'Asp_DL'])),]
maxes <- rbind(maxes,dataOne[(dataOne[,'Glu_DL']==max(dataOne[,'Glu_DL'])),])
maxes <- rbind(maxes,dataOne[(dataOne[,'Ala_DL']==max(dataOne[,'Ala_DL'])),])
maxes <- rbind(maxes,dataOne[(dataOne[,'Val_DL']==max(dataOne[,'Val_DL'])),])
maxes <- rbind(maxes,dataOne[(dataOne[,'Phe_DL']==max(dataOne[,'Phe_DL'])),])


#########################################################################################
##	AMINO ACID D/L RESIDUALS LIST
#########################################################################################

quants <- matrix(nrow=nrow(dataToUse),ncol=5)
for (i in 1:nrow(resid)) {
	quants[i,] <- quantile((resid[i,]), c(0.025,0.25,0.5,0.75,0.975), na.rm=TRUE)
}

par(mfcol=c(1,1))
#iniA4 <- iniA3[1:3]
iniA4 <- iniA3
yInt <- 0.005

plot(c(1,nrow(resid)),c(min(resid, na.rm=TRUE)-yInt*(length(iniA4)+1),max(resid, na.rm=TRUE)), type='n', xlab='',xaxt='n')
#plot(c(1,nrow(resid)),c(-0.08,0.08), type='n', xlab='',xaxt='n')
axis(1,at=seq(1,nrow(resid), by=1),labels=dataToUse[,'specimen_no'], las=3)
lines(c(-1,nrow(resid))+2,c(0,0))
segments(x0=seq(1,nrow(resid), by=1),y0=quants[,1],y1=quants[,5], col='grey')
segments(x0=seq(1,nrow(resid), by=1),y0=quants[,2],y1=quants[,4], lwd=5, col='grey')
for (x in iniA4) {
	for (y in iniA4) {
		points(seq(1,nrow(resid), by=1),(resid[,paste(x,y,sep='')]), pch=19, cex=0.5)
		text(seq(1,nrow(resid), by=1),(resid[,paste(x,y,sep='')]), paste(x,y,sep=''), cex=0.5, pos=4)
}}
#iniA3 <- c('Asp','Glu','Ala','Val','Phe')

yMin <- min(resid, na.rm=TRUE)

for (x in iniA4) {
	yMin <- yMin - yInt
	text(-1,yMin,x, cex=0.6)
	text(seq(1,nrow(resid), by=1), yMin, round(dataToUse[,paste(x,'DL',sep='_')],3), cex=0.3)
}
yMin <- yMin - yInt
text(-1,yMin,'14C', cex=0.6)
text(seq(1,nrow(resid), by=1), yMin,dataToUse[,agesToUse[2]], cex=0.3)

text(0, max(resid, na.rm=TRUE), taxon, cex=2, pos=4)


##	AMINO ACID D/L DETAILED CORRELATIONS
#########################################################################################

par(mfcol=c(length(iniA3),length(iniA3)), cex=0.5, mar=c(1,1,1,1))
par(mfcol=c(length(iniA4),length(iniA4)), cex=0.5, mar=c(1,1,1,1))

resid <- matrix(nrow=nrow(dataOne),ncol=length(iniA3)*length(iniA3))
colname <- list()
i <- 0
for (x in iniA3) {
	for (y in iniA3) {
		i<-i+1
		colname[i] <- paste(x,y,sep='')
	}
}
colnames(resid) <- colname

ChooseMe <- c(10098,10111,10114,10119,10120,10123,10132,10143,10139)
ChooseMe <- c(10103,10104,10109,10114,10116,10136,10138,10120,10094,10097,10111,10139,10140)
ChooseMe <- c(ChooseMe,paste(ChooseMe,'R',sep=''))
ChooseMe <- c(ChooseMe,'10126A','10126B','10126C','10126D')
# NO 10122,10123,10106

#for (t in c(50,84)) {
#	tData <- dataOne[(dataOne[,'top']==50),]
	tData <- rbind(dataOne[(dataOne[,'top']==50),], dataOne[(dataOne[,'top']==84),])
	tData2 <- subset(tData, Glu_DL < 0.12 & Glu_DL>0.08)
	ChooseMe2 <- c(10328,10335,10321)
	ChooseMe2 <- c(10317,10335,10321,10320,10319,'10321D','10321E')

ChooseMe3 <- c(10094,10104,10111,10114,10136,10138,10139)
ChooseMe3 <- c(ChooseMe3,paste(ChooseMe3,'R',sep=''))


for (x in iniA4) {
	for (y in iniA4) {
		
		if (x==y) {
			plot(1:1, type='n', xaxt='n', yaxt='n', ylab='', xlab='')
			text(1,1.1,paste('[','DL ',x,']', sep=''), cex=3)
			text(1,0.9,taxon, cex=2, font=4)
			text(1,0.75,'Watson Bay Core', cex=1.3)
		} else {
			
			xCol <- paste(x,'DL',sep='_')
			yCol <- paste(y,'DL',sep='_')
			
			yLim <- c(min(dataOne[,yCol]),max(dataOne[,yCol])*1.1)
			xLim <- c(min(dataOne[,xCol]),max(dataOne[,xCol])*1.1)
			
			plot(dataOne[,xCol], dataOne[,yCol], xlab='', ylab='', xlim=xLim, ylim=yLim,pch=21, col='white', bg='grey', lwd=0.2, cex=0.5)
		
			a <- cor.test(dataOne[,xCol], dataOne[,yCol], method='spearman', exact=FALSE)
			b <- lm(dataOne[,yCol] ~ dataOne[,xCol])
			resid[, paste(x,y,sep='')] <- b$residuals
			abline(b)		
		
			text(max(dataOne[,xCol])*1.0, min(dataOne[,yCol])*1.2, paste('rho=',round(a$estimate,2)), pos=2)
			text(min(dataOne[,xCol])*1.2, max(dataOne[,yCol])*0.9, paste('r2=',round(summary(b)$r.squared,2)), pos=4)

			points(tData[,xCol], tData[,yCol],pch=21,col='orange', bg='orange', cex=0.7)
#			text(tData2[,xCol], tData2[,yCol], tData2[,'ual_num'], cex=0.7, pos=4)
			for (sp in ChooseMe2) {
				points(dataOne[(dataOne[,'ual_num']==sp),xCol], dataOne[(dataOne[,'ual_num']==sp),yCol],pch=21,col='orange', bg='black', cex=0.7)
				text(dataOne[(dataOne[,'ual_num']==sp),xCol], dataOne[(dataOne[,'ual_num']==sp),yCol], dataOne[(dataOne[,'ual_num']==sp),'ual_num'], cex=0.7, pos=4)				
			}			
			for (sp in ChooseMe) {
				points(dataOne[(dataOne[,'ual_num']==sp),xCol], dataOne[(dataOne[,'ual_num']==sp),yCol],pch=21,col='blue', bg='black', cex=0.7)
				text(dataOne[(dataOne[,'ual_num']==sp),xCol], dataOne[(dataOne[,'ual_num']==sp),yCol], dataOne[(dataOne[,'ual_num']==sp),'ual_num'], cex=0.7, pos=4)				
			}
			for (sp in ChooseMe3) {
				points(dataOne[(dataOne[,'ual_num']==sp),xCol], dataOne[(dataOne[,'ual_num']==sp),yCol],pch=21,col='red', bg='red', cex=1.0)
			}
		}
	}
}



##	LIST
#########################################################################################

par(mfcol=c(1,1))
iniA4 <- iniA3[1:3]

plot(c(1,nrow(resid)),c(min(resid,na.rm=TRUE),max(resid, na.rm=TRUE)), type='n')
plot(c(1,nrow(resid)),c(-0.08,0.08), type='n', xlab='',xaxt='n')
axis(1,at=seq(1,nrow(resid), by=1),labels=dataOne[,'ual_num'], las=3)
lines(c(-1,nrow(resid))+2,c(0,0))
segments(x0=seq(1,nrow(resid), by=1),y0=quants[,1],y1=quants[,5], col='grey')
segments(x0=seq(1,nrow(resid), by=1),y0=quants[,2],y1=quants[,4], lwd=5, col='grey')
for (x in iniA4) {
	for (y in iniA4) {
		points(seq(1,nrow(resid), by=1),(resid[,paste(x,y,sep='')]), pch=19, cex=0.5)
		text(seq(1,nrow(resid), by=1),(resid[,paste(x,y,sep='')]), paste(x,y,sep=''), cex=0.5, pos=4)
}}

quants <- matrix(nrow=nrow(dataOne),ncol=5)

for (i in 1:nrow(resid)) {
	quants[i,] <- quantile((resid[i,]), c(0.025,0.25,0.5,0.75,0.975), na.rm=TRUE)
}

ChooseMe <- c(10098,10111,10114,10119,10120,10123,10132,10143,10139)
ChooseMe <- c(10103,10104,10109,10114,10116,10136)









DL <- 'Ala_DL'
xLab <- expression('uncalibrated Ala D/L'^{2})
plot(sp2data[,DL]^2,sp2data[,'top'], ylim=c(170,0), ylab='uncorrected core depth (m)', xlab=xLab, main='Watsons Bay (AUSCDT) Core', type='n', yaxt='n')
axis(2,seq(0,180,by=30), seq(0,1.8,by=0.3), las=1)

for (d in TOPS) {
	a <- quantile(sp2data[(sp2data[,'top']==d),DL]^2, c(0.025,0.25,0.5,0.75,0.975), rm.na=TRUE)
#	lines(c(a[1],a[5]),c(d,d))
	bot <- as.numeric(d)-2
	top <- as.numeric(d)+2
	polygon(c(a[2],a[4],a[4],a[2]),c(bot,bot,top,top), col='lightgrey')
#	points(a[3],d,pch=22, col='blue', bg='blue')
	lines(c(a[3],a[3]),c(top+2,bot-2), lwd=2)
}
points(sp2data[,DL]^2,sp2data[,'top'], pch=21,col='white',bg='black')
text(min(sp2data[,DL]^2),max(sp2data[,'top'])+10,'Fulvia tenuicostata',font=4, pos=4)


DL <- 'Ala_DL'
xLab <- expression('uncalibrated Ala D/L'^{2})
plot(sp1data[,DL]^2,sp1data[,'top'], ylim=c(170,0), ylab='uncorrected core depth (m)', xlab=xLab, main='Watsons Bay (AUSCDT) Core', type='n', yaxt='n')
axis(2,seq(0,180,by=30), seq(0,1.8,by=0.3), las=1)

for (d in TOPS) {
	if (length(sp1data[(sp1data[,'top']==d),DL])>0){
	a <- quantile(sp1data[(sp1data[,'top']==d),DL]^2, c(0.025,0.25,0.5,0.75,0.975), na.rm=TRUE)
#	lines(c(a[1],a[5]),c(d,d))
	bot <- as.numeric(d)-2
	top <- as.numeric(d)+2
	polygon(c(a[2],a[4],a[4],a[2]),c(bot,bot,top,top), col='lightgrey')
#	points(a[3],d,pch=22, col='blue', bg='blue')
	lines(c(a[3],a[3]),c(top+2,bot-2), lwd=2)
}}
points(sp1data[,DL]^2,sp1data[,'top'], pch=21,col='white',bg='black')
text(min(sp1data[,DL]^2),max(sp1data[,'top'])+10,'Timoclea (Chioneryx) cardioides',font=4, pos=4)

par(mfcol=c(1,3))

DL <- 'Ala_DL'
for (DL in c('Asp_DL','Glu_DL','Ala_DL')) {
xLab <- expression('uncalibrated Ala D/L'^{2})
plot(sp3data[,DL]^2,sp3data[,'top'], ylim=c(170,0), ylab='uncorrected core depth (m)', xlab=xLab, main='Watsons Bay (AUSCDT) Core', type='n', yaxt='n')
axis(2,seq(0,180,by=30), seq(0,1.8,by=0.3), las=1)

for (d in TOPS) {
	if (length(sp1data[(sp3data[,'top']==d),DL])>0){
	a <- quantile(sp3data[(sp3data[,'top']==d),DL]^2, c(0.025,0.25,0.5,0.75,0.975), na.rm=TRUE)
#	lines(c(a[1],a[5]),c(d,d))
	bot <- as.numeric(d)-2
	top <- as.numeric(d)+2
	polygon(c(a[2],a[4],a[4],a[2]),c(bot,bot,top,top), col='lightgrey')
#	points(a[3],d,pch=22, col='blue', bg='blue')
	lines(c(a[3],a[3]),c(top+2,bot-2), lwd=2)
	}
}
points(sp3data[,DL]^2,sp3data[,'top'], pch=21,col='white',bg='black')
text(min(sp3data[,DL]^2),max(sp3data[,'top'])+10,'Peronella peronii',font=4, pos=4)
}
	