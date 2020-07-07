# BIN NAMES / DATES
binNames <- c('Cam1','Cam2','Cam3','Cam4','Ord1','Ord2','Ord3','Ord4','Ord5','Sil1','Sil2','Dev1','Dev2','Dev3','Dev4','Dev5','Car1','Car2','Car3','Car4','Car5','Per1','Per2','Per3','Per4','Tri1','Tri2','Tri3','Tri4','Jur1','Jur2','Jur3','Jur4','Jur5','Jur6','Cre1','Cre2','Cre3','Cre4','Cre5','Cre6','Cre7','Cre8','Cen1','Cen2','Cen3','Cen4','Cen5','Cen6')

binsR <-c(0, 11.6, 23.0, 33.9, 40.4, 55.8, 65.5, 70.6, 83.5, 93.5, 99.6, 112.0, 125.0, 136.4, 145.5, 150.8, 164.7, 171.6, 183.0, 189.6, 199.6, 216.5, 228.0, 245.0, 251.0, 260.4, 270.6, 284.4, 299.0, 306.5, 318.1, 336.0, 349.5, 360.7, 367.1, 383.7, 391.9, 409.1, 418.1, 428.2, 443.7, 449.5, 460.5, 466.0, 479.0, 490.0, 501.0, 513.0, 530.0,542.0)

bins <- vector(length=length(binsR))

for (i in 1:(length(bins))) {
	bins[i] <- binsR[(length(bins)-i+1)]
}

# CALCULATE BIN MID POINTS
binMeta <- matrix(nrow=length(bins)-1,ncol=4)
colnames(binMeta) <- c('name','boundary','length','mids')

binMids <- vector(mode='numeric', length=length(bins))

for (i in 1:(length(bins)-1)) {
	binMids[i] <- bins[i] + ((bins[i+1] - bins[i])/2)
	
	binMeta[i,'mids'] <- binMids[i]
	binMeta[i,'length'] <- bins[(i+1)] - bins[i]
	binMeta[i,'boundary'] <- bins[i+1]
	binMeta[i,'name'] <- binNames[length(binNames)-i+1]
	
}

binMids <- sort(binMids, decreasing=TRUE)

# PERIOD BREAKS
peri <-c(0,65.5,145.5,199.6,251.0,299.0,360.7,418.1,443.7,490.0,542.0)
peri_n <-c('T','K','J','Tr','P','C','D','S','O','Cm')

# ERA BREAKS
big <-c(0,65.5,251.0,542.0)
big_n <-c('Cz','Mz','Pz')
