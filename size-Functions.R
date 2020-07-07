#########################################################################################
plotSizeCrossCorrelations <- function(areaRaw, fileName) {
	
	if (!is.na(fileName))
		pdf(paste(PATH.OUT,fileName,sep='/'), width=pageWidthTwo, height=pageWidthTwo)
#	par(mfcol=c(2,2), cex=0.5, mar=c(1,1,1,1), mgp=c(1,0.5,0))
	par(mfcol=c(2,2))
	
	plot(areaRaw[,'x_mm'],areaRaw[,'y_mm'])
	plot(areaRaw[,'x_mm'],areaRaw[,'z_mm'])
	plot(areaRaw[,'y_mm'],areaRaw[,'z_mm'])
	plot(areaRaw[,'mass_mg'],areaRaw[,'gm_size'], log='xy')

	if (!is.na(fileName))
		dev.off()
	
}


#########################################################################################
plotSizeMassCorrelations <- function(areaRaw, fileName) {
	
	if (!is.na(fileName))
		pdf(paste(PATH.OUT,fileName,sep='/'), width=pageWidthTwo, height=pageWidthTwo)
#	par(mfcol=c(2,2), cex=0.5, mar=c(1,1,1,1), mgp=c(1,0.5,0))
	par(mfcol=c(2,2))
	
	a <- lm(log(areaRaw[,'mass_mg']) ~ log(areaRaw[,'gm_size']))
	plot(a)

	if (!is.na(fileName))
		dev.off()
	
	return (a)
}


#########################################################################################
plotSizeBySample <- function(areaRaw, colGrp, fileName) {

	GROUPS <- unique(as.factor(areaRaw[,colGrp]))

	if (!is.na(fileName))
		pdf(paste(PATH.OUT,fileName,sep='/'), width=pageWidthTwo, height=pageWidthTwo)
	
	par(mfcol=c(2,1), cex=0.5, mar=c(2,2,1,1), mgp=c(1,0.5,0))
		
	for (j in 1:2) {
		
		if (j ==1) {
			size <- 'gm_size'
			yLab <- 'size (mm)'
		}
		if (j == 2) {
			size <- 'mass_mg'
			yLab <- 'mass (mg)'	
		}
		
		areaPlot <- areaRaw[!is.na(areaRaw[,size]) & (areaRaw[,size] > 0),]
		
		plot(as.factor(areaPlot[,colGrp]), areaPlot[,size], log='y', ylab=yLab)
		points(as.factor(areaPlot[,colGrp]), areaPlot[,size], pch=19)
		
		sq <- quantile(areaPlot[,size], na.rm=TRUE, qSD)
		for (i in 1:length(sq))
			lines(c(0,length(GROUPS)+1),c(sq[i],sq[i]), col='grey')
		
	}

	if (!is.na(fileName))
		dev.off()
}

#########################################################################################
sumByGroup <- function (dataOne, sumCols, groupBy) {

	GROUPS <- unique(dataOne[,groupBy])

	for (g in GROUPS) {

		specRows <- dataOne[(dataOne[,groupBy]==g),]

		if (nrow(specRows)>1) {
			for (w in sumCols) {
				specRows[1,w] <- sum(specRows[,w],na.rm=TRUE)
			}	
		}
		if (g == GROUPS[1])
			dataSpec <- dataOne[1,]
		else 
			dataSpec <- rbind(dataSpec,specRows[1,])
	}

	return(dataSpec)
}

#########################################################################################
specimenCollectionSize <- function(dataOne) {
	
	specs <- unique(dataOne[,'specimen_no']) 
	
	for (s in specs) {
		specRows <- dataOne[(dataOne[,'specimen_no'] == s) & !is.na(dataOne[,'x_mm']),]

		if (nrow(specRows)>1) {
			firstMeas <- min(as.Date(specRows[,'date_measured'],'%Y-%m-%d'),na.rm=TRUE)
			firstData <- specRows[(as.Date(specRows[,'date_measured'],'%Y-%m-%d') == firstMeas),]
		} else {
			firstData <- specRows
		}
		
		if (s == specs[1])
			uniqueData <- firstData
		else
			uniqueData <- rbind(uniqueData, firstData)
		
	}
	
	return(uniqueData)
}

#########################################################################################
estimateMassFromSize <- function(sourceData, sinkData) {

	taxaSource <- unique(sourceData[,'taxon_no']) 
	taxaSink <- unique(sourceData[,'taxon_no']) 
	taxa <- union(taxaSource, taxaSink)

	##	MAKE SUMMARY MATRIX FOR TAXON SIZE RELATIONS
	taxCol <- c('taxon_no','nSpec','r2','m','b','mult')
	taxData <- matrix(ncol=length(taxCol),nrow=length(taxa))
	colnames(taxData) <- taxCol
	
	for (t in 1:length(taxa)) {
		
		##	GET REQUIRED DATA SUBSET
		tData <- sourceData[(sourceData[,'taxon_no']== taxa[t]),]
		tData <- tData[(!is.na(tData[,'x_mm']) & !is.na(tData[,'y_mm']) & !is.na(tData[,'z_mm']) & !is.na(tData[,'mass_mg'])),]	
		
		##	TAXON META DATA
		taxData[t,'taxon_no'] <- taxa[t]
		taxData[t,'nSpec'] <- nrow(tData)
		
		##	ONLY FIT N > 5 
		if (nrow(tData) > 5) {
						
			## 	ZERO X ZERO IS THE CORRECT INTERCEPT
			a <- lm (I(mass_mg) ~ I(gm_size^3) + 0, data=tData)
			taxData[t,'r2'] <- summary(a)$r.squared
			taxData[t,'m'] <- summary(a)$coefficients[1,1]
			taxData[t,'b'] <- 0
	
			##	NEED TO ADJUST FOR LIVE SPECIMENS HAVING TWO VALVES
			if ("right valve" %in% unique(tData[,'measured_object'])) {
				taxData[t,'mult'] <- 2
			} else {
				taxData[t,'mult'] <- 1
			}
		}
	}
	
	for (t in 1:length(taxa)) {

		tData <- sinkData[(sinkData[,'taxon_no']== taxa[t]),]
		
		if (taxData[t,'r2'] > 0.9)			
			tData[is.na(tData[,'mass_mg']),'mass_mg'] <- (tData[is.na(tData[,'mass_mg']),'gm_size']^3) * taxData[t,'m']
		
		if (t == 1) {
			bData <- tData
		} else {
			bData <- rbind(bData,tData)
		}
	}

	return(bData)
}

#########################################################################################
massSummaryMatrix <- function(rawData, pKey, keyCols, rootName, fracs, upKeys) {
	
	if (missing(upKeys)) {
		upKeys <- sort(unique(rawData[, pKey]))
	}
	
	if (is.na(keyCols[1])) {
		sumCol <- c(pKey, paste(rootName,fracs,sep='.'))
	} else {
		sumCol <- c(keyCols, pKey, paste(rootName,fracs,sep='.'))
	}
	sumData <- matrix(ncol=length(sumCol),nrow=length(upKeys))
	colnames(sumData) <- sumCol
	
	sumData[,pKey] <- upKeys
	
	for (u in upKeys) {
		
		kData <- rawData[(rawData[,pKey] == u) & !is.na(rawData[,pKey]),]
		if (nrow(kData) > 0) {
			if ('site_no' %in% c(pKey,keyCols))
				sumData[(sumData[,pKey] == u),'site_no'] <- kData[1,'site_no']
			if ('sample_no' %in% c(pKey,keyCols))
				sumData[(sumData[,pKey] == u),'sample_no'] <- kData[1,'sample_no']
			if ('layer_no' %in% c(pKey,keyCols))
				sumData[(sumData[,pKey] == u),'layer_no'] <- kData[1,'layer_no']
		
			# SUM AND COVERT MG TO GRAMS
			for (f in fracs) {
				sumData[(sumData[,pKey] == u),paste(rootName,f,sep='.')] <- sum(kData[(kData[,'sieve_size'] == f),'mass_mg'], na.rm=TRUE)/1000
			}
		}
	}
	
	return(sumData)
}

#########################################################################################
combineSummaryMatrices <- function(matrix1, matrix2,pKey) {
	
	orderMatch <- TRUE
	
	if (nrow(matrix1) == nrow(matrix2)) {
		for (i in 1:nrow(matrix1)) {
			if ((matrix1[i,pKey] != matrix1[i,pKey]) & (orderMatch))
				orderMatch <- FALSE	
		}
		if (orderMatch)
			return(cbind(matrix1,matrix2))
	}
	return(FALSE)
}


#########################################################################################
addProportionColumns <- function(matrix, taxon, fracs) {
	
	colHead <- colnames(matrix)
	nHead <- ''
	
	if ( 16 %in% fracs) {
		p16 <- matrix[,paste(taxon,'16',sep='.')]/matrix[,'mass.16']
		matrix <- cbind(matrix,p16)
		nHead <- c(nHead, paste('p',taxon,'16',sep='.'))
	}

	if ( 8 %in% fracs) {
		p8 <- (matrix[,paste(taxon,'16',sep='.')]+matrix[,paste(taxon,'8',sep='.')]) / (matrix[,'mass.16']+matrix[,'mass.8'])
		matrix <- cbind(matrix,p8)
		nHead <- c(nHead, paste('p',taxon,'8',sep='.'))
	}

	if ( 4 %in% fracs) {
		p4 <- (matrix[,paste(taxon,'16',sep='.')]+matrix[,paste(taxon,'8',sep='.')]+matrix[,paste(taxon,'4',sep='.')]) / (matrix[,'mass.16']+matrix[,'mass.8']+matrix[,'mass.4'])
		matrix <- cbind(matrix,p4)
		nHead <- c(nHead, paste('p',taxon,'4',sep='.'))
	}
	nHead <- nHead[-1]
	colnames(matrix) <- c(colHead,nHead)

	return(matrix)
}
