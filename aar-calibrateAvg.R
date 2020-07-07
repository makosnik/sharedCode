###################################################
## NEEDS: datesToUse
## RETURNS: datesCalibrated / SAVES: 'taxon'-CalibratedData.csv
###################################################


# NOTIFY SCREEN
###################################################
print('start aar-calibrateAvg.R')

# VARIABLES ASSUMED TO BE SET
###################################################
# LIVING <- "include"
# taxon

c14Cage <- 'X2sMedian'
c14Clwr <- 
c14Cupr <- 

# LOAD EXTERNAL FILES
###################################################
source('../sharedCode/aar-calibrationFunctions.R')
source('../sharedCode/regress.R')
source('../sharedCode/histogram.R')


# SET - IMPORTANT CONSTANTS
###################################################
FILEROOT <- paste(FILEOUT,taxon,sep='')

DATEZERO <- 1950
dateMod <- 56

# LIVING DATA Options
##############################
if (LIVING == 'exclude') {

	print ('excluding - live data')
	datesToUse[((datesToUse[, c14Cage] < -50) & (!is.na(datesToUse[, c14Cage]))),c(c14Cage, c14Clwr, c14Cupr)]<- NA
	FILEROOT <- paste(FILEROOT,"LiveExcluded",sep='-')

} else if (LIVING == 'origin') {

	print ('using - live data to force origin')
	FILEROOT <- paste(FILEROOT,"LiveOrigin",sep='-')

} else {

	FILEROOT <- paste(FILEROOT,"LiveIncluded",sep='-')

}

###################################################
## START - CALIBRATE AMINO ACID DATES USING CARBON
###################################################

tmp <- ''
#print(paste('DATEZERO =',DATEZERO))

# START FORCE LIVING ORIGIN Hack
##############################
if (LIVING == 'origin') {

	# CHANGE FROM AD 1950 TO AD XXXX
	if (dateMod == 0)
		dateMod <- min(datesToUse[,'X2sMedian'], na.rm=TRUE)
	if (dateMod < 0)
		dateMod <- abs(dateMod)
	else {
		dateMod <- 56
#			print (paste("XXXXXXXXX USING ASSUMED dateMod XXXXXXXXX",dateMod))
	}

	# REMOVE LIVING SHELLS ABSOLUTE AGES
	datesToUse[((datesToUse[,'X2sMedian'] < -50) & (!is.na(datesToUse[,'X2sMedian']))),c('X2sMedian','X2sYng','X2sOld')]<- NA

	# SET AD XXXX TO 0 NOT AD 1950
	datesToUse[,'X2sMedian'] <- datesToUse[,'X2sMedian'] + dateMod
	datesToUse[,'X2sYng'] <- datesToUse[,'X2sYng'] + dateMod
	datesToUse[,'X2sOld'] <- datesToUse[,'X2sOld'] + dateMod

	DATEZERO <- 1950 + dateMod
	print (paste('using - date relative to AD',DATEZERO))
}
# FINISH FORCE LIVING ORIGIN Hack

# DATA - ONLY SPECIMENS WITH CARBON DATES	
c14 <- datesToUse[!is.na(datesToUse[,'X2sMedian']),]



# FILE - OPEN CALIBRATION PLOT FILE
##############################
postscript(paste(FILEROOT,"Plot.eps",sep="-"))
par(mar=c(4,4,2,2), las=0, mgp=c(2,1,0))
plotSize <- 7
rows <- 3
cols <- 2
nplots <- rows*cols
layout(matrix(seq(nplots), rows, cols, byrow=FALSE), width=rep(lcm(plotSize),cols), height=rep(lcm(plotSize),rows), respect=TRUE)
plotLetter <- 0


# PROCEED ONLY IF THE TAXON HAS MORE THAN FOUR DATES
##############################
if (nrow(c14) > 4) {

	# EACH AMINO ACID HAS OWN CALIBRATION CURVE
	##############################
	for (amino in 1:length(colR)) 
		if (!is.na(mean(datesToUse[,colR[amino]],rm.na=true))) {
	
		print(paste('start calibrating',taxon,'using', namA[amino],' [aarCalibrateAvg.R]'))
	
		
		if (LIVING == 'origin') {
	
			# GET LIVING D/L VALUES
			# livingDL <- datesToUse[((datesToUse[,'X2sMedian'] < -50) & (!is.na(datesToUse[,'X2sMedian']))),c(colR)]
		#	print ('using - live D/L as live value')
		
			# SUBTRACT LIVING D/L FROM ALL SPECIMENS
			datesToUse[,colR[amino]] <- datesToUse[,colR[amino]] - min(datesToUse[,colR[amino]])
			print ('using - minimmum D/L as live value')
	
		}
			
		# PERFORM WEIGHTED CALIBRATIONS...
		##############################
		DATA <- ''
		
		ageMin <- -100
#		ageMax <- max(datesToUse[,'X2sOld'],na.rm=TRUE)
		ageMax <- NA
#		if (taxon =='Tellinid')
#			ageMax <- 2500
#		if (taxon == 'Liloa')
#			ageMax <- 3200
		
		w1 <- aaCalCurveOptim(c(0.4,6), ageMin, ageMax, data14C=datesToUse[,c('X2sMedian','X2sYng','X2sOld')], dataAmino=datesToUse[,colR[amino]], nameAmino=iniA[amino], nameSpecimen=taxon, DATEZERO=DATEZERO, LIVING=LIVING)
		plotLetter <- plotLetter + 1
		title (main=paste(LETTERS[plotLetter],'.'),adj=0)
	
		title (main=namA[amino], adj=0.5, font.main=1)
		
		# DATA - MAKE FINAL DATA
		##############################
	#	if (LIVING == 'origin') {
	#		w1 <- w1 - dateMod
	#	}
	
		# INFER AGES TO THE NEAREST YEAR - NOT MORE PRECISION THAN THAT!
		w1 <- round(w1,0)
	
		# DATA - GIVE CALIBRATED DATE COLUMNS MEANINGFUL NAMES
		colnames(w1) <- c(paste(iniA[amino],'AgeFit',sep='_'),paste(iniA[amino],'AgeLwr',sep='_'),paste(iniA[amino],'AgeUpr',sep='_'))
	
		datesToUse <- cbind(datesToUse,w1)
	}
	
	
	ageCols <- vector(mode='logical', length(colR))
	for (amino in 1:length(colR))
		ageCols[amino] <- paste(iniA[amino],'AgeFit',sep='_')
	
	ageCols <- intersect(ageCols,colnames(datesToUse))
	
	ageUncert <- matrix(nrow=nrow(w1),ncol=length(colR))
	mn <- vector(mode='logical', length(colR))
	#ageCols <- vector(mode='logical', length(colR))
	for (amino in 1:length(colR))
		mn[amino] <- paste(iniA[amino],'AgeUnc',sep='_')
	colnames(ageUncert)<-mn
	
	big <- 10000
	
	for (amino in 1:length(colR))
	if (!is.na(mean(datesToUse[,colR[amino]],rm.na=true))){
		ageCols[amino] <- paste(iniA[amino],'AgeFit',sep='_')
		ageUncert[,paste(iniA[amino],'AgeUnc',sep='_')] <- (datesToUse[,paste(iniA[amino],'AgeUpr',sep='_')]+big) - (datesToUse[,paste(iniA[amino],'AgeLwr',sep='_')]+big)
	}

}
#datesCalibrated <- cbind(datesToUse, weightedAgeMean, weightedAgeSdev)
datesCalibrated <- datesToUse

# FILE - CLOSE CALIBRATION PLOT FILE
##############################
dev.off()
# DATA - SAVE CALIBRATED DATA
###################################################
# CLEANEDBY, 
save(datesCalibrated, colR, iniA, taxon, LIVING, file=paste(FILEDAT,taxon,"-aarCalibratedAvg.Rdata",sep=''))

# FILE - DUMP REJECT SEQUENCE AS LIST
of <- paste(FILEROOT,"CalibratedData.csv",sep='-')
write.csv(datesCalibrated, file=of)


###################################################
## END - CALIBRATE AMINO ACID DATES USING CARBON
###################################################

print('finish aar-calibrateAvg.R')
