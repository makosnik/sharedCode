
##	ADD MASS / AREA COLUMNS TO SUMMARY TABLE
#######################################################
addMassAreaColumns <- function(dataMatrix,colMassRoot,fractions) {

	dataMatrixCols <- colnames(dataMatrix)

	for (f in fractions) {
		massArea <- dataMatrix[,paste(colMassRoot,f,sep='.')] / dataMatrix[,'sample_size']
		dataMatrix <- cbind(dataMatrix, massArea)
	}
	
	colnames(dataMatrix) <- c(dataMatrixCols,paste(colMassRoot,'Area',fractions,sep='.'))
	
	return(dataMatrix)
}

##	ADD TIME-AVERAGING COLUMNS TO SUMMARY TABLE
#######################################################
addTimeAveragingColumns <- function(dataMatrix) {

	ta1 <- dataMatrix[,paste('age', qSD1[3],sep='.')] - dataMatrix[,paste('age', qSD1[1],sep='.')]
	ta2 <- dataMatrix[,paste('age', qSD2[3],sep='.')] - dataMatrix[,paste('age', qSD2[1],sep='.')]
	taIQ <- dataMatrix[,paste('age', qIQ0[3],sep='.')] - dataMatrix[,paste('age', qIQ0[1],sep='.')]

	dataMatrix <- cbind(dataMatrix,ta1,ta2,taIQ)
	
	return(dataMatrix)
}


##	ADD CARBONATE PRODUCTION COLUMNS TO SUMMARY TABLE
#######################################################
addCarbonateProductionColumns <- function(dataMatrix, colMARoot, pmToUse, taToUse, sOut, fracs) {

	dataMatrixCols <- colnames(dataMatrix)

	colsMA <- paste(colMARoot, fracs, sep='.')

	for (m in colsMA) {
		carbProd <- ((dataMatrix[,m]) * pmToUse) / dataMatrix[,taToUse]
		dataMatrix <- cbind(dataMatrix, carbProd)	
	}

	colnames(dataMatrix) <- c(dataMatrixCols,paste(sOut,'CP',fracs,taToUse,sep='.'))
	
	return(dataMatrix)
}


#######################################################
createAgesSummaryTable <- function(rawData, pKey, quants) {
	
	pKeys <- sort(unique(rawData[, pKey]))
	
	stdCols <- c('top','depth','sample_size')
	
	sColNames <- c(pKey,stdCols,'nAges',paste('age',quants,sep='.'))
	sMatrix <- matrix(nrow=length(pKeys),ncol=length(sColNames))
	colnames(sMatrix) <- sColNames

	sMatrix[,pKey] <- pKeys
	
	for (k in 1:length(pKeys)) {
		kData <-rawData[(rawData[,pKey] == pKeys[k]),]

		## ROW SPECIFIC DATA
		for (col in stdCols) {
			if (col %in% colnames(rawData))
				sMatrix[(sMatrix[,pKey] == pKeys[k]),col] <- kData[1,col]		
		}
		
		## AGES
		ages <- kData[, agesToUse[2]]
		ages <- ages[!is.na(ages)]
		a <- quantile(ages, quants)
		sMatrix[k,'nAges'] <- length(ages)
		sMatrix[k,paste('age', quants, sep='.')] <- a
		
	}
	return(sMatrix)
}

#######################################################
createSedimentSummaryTable <- function(rawData, pKey, fracs) {
	
	pKeys <- sort(unique(rawData[, pKey]))
	
	stdCols <- c('top','depth')
	
	sColNames <- c(pKey,stdCols,'nKey','sample_size',paste('mass', fracs,sep='.'))
	sMatrix <- matrix(nrow=length(pKeys),ncol=length(sColNames))
	colnames(sMatrix) <- sColNames

	sMatrix[,pKey] <- pKeys
	
	for (k in 1:length(pKeys)) {
		kData <-rawData[(rawData[,pKey] == pKeys[k]),]

		## ROW SPECIFIC DATA
		for (col in stdCols) {
			if (col %in% colnames(rawData))
				sMatrix[(sMatrix[,pKey] == pKeys[k]),col] <- kData[1,col]		
		}

		sMatrix[k,'nKey'] <- length(unique(kData[, pKey]))
		
		##	FOR EACH FRACION COMBINE SEDIMENT DATA INTO SITE SUMMARY
		for (m in fracs) {
			
			siteSize <- kData[(kData[,'sieve_size']==m),]
	
			if (nrow(siteSize) >= sMatrix[k,'nKey']) {
				if (is.na(sMatrix[k,'sample_size'])) {
					sMatrix[k,'sample_size'] <- sum(siteSize[,'sample_size'])
				} else {
					if (sMatrix[k,'sample_size'] != sum(siteSize[,'sample_size']))
						print('error')
				}	
				sMatrix[k,paste('mass',m,sep='.')] <- sum(siteSize[,'sediment_mass_g'], na.rm=TRUE) - sum(siteSize[,'noncarbonate_mass_g'], na.rm=TRUE)
			}
		}
		
	}
	return(sMatrix)
}

#######################################################################################################
combineSummaryTables <- function(matrix1, matrix2, pKey) {
	
	mk1 <- matrix1[,pKey]
	mk2 <- matrix2[,pKey]
	uKeys <- intersect(mk1, mk2)
	
	mh1 <- colnames(matrix1)
	mh2 <- colnames(matrix2)
	mhi <- unique(c(mh1, mh2))
	
	cMatrix <- matrix(nrow=length(uKeys),ncol=length(mhi))
	colnames(cMatrix) <- mhi
	cMatrix[,pKey] <- uKeys
	
	for (k in uKeys) {
		cMatrix[(cMatrix[,pKey] == k), mh1] <- matrix1[(matrix1[,pKey] == k),]
		cMatrix[(cMatrix[,pKey] == k), mh2] <- matrix2[(matrix2[,pKey] == k),]
	}
	
	return(cMatrix)
}


#######################################################################################################
addTaxonMassSummaryColumnsByLayer <- function(tax, taxData, keyName, keyValues, fracs) {

	sColNames <- c(keyName, paste('mass',fracs,tax,sep='.'))
	sMatrix <- matrix(nrow=length(keyValues),ncol=length(sColNames))
	colnames(sMatrix) <- sColNames

	for (l in 1:length(keyValues)) {
		sMatrix[l, keyName] <- keyValues[l]
		for (f in fracs) {
			colName <- paste('mass',f,tax,sep='.')
			taxFrac <- taxData[(taxData[,keyName]==keyValues[l]) & (taxData[,'sieve_size']==f),]
			sMatrix[l, colName] <- sum(taxFrac[,'mass_mg'])/1000
		}
	}
	return(sMatrix)
}



