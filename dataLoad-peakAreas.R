###################################################
# DATA - FILE PREPARATION
###################################################
# SQL TO GET DATA FOR INPUT FILE
# SELECT mak_taxonomy.names.name, aarPeakAreas.*, c14Dates.*, measure2.*, layers.* FROM aarPeakAreas LEFT JOIN c14Dates USING (specimen_no) LEFT JOIN measure2 USING (specimen_no) LEFT JOIN specimens USING (specimen_no) JOIN fractions USING (fraction_no) JOIN layers USING (layer_no) JOIN mak_taxonomy.names USING (taxon_no) INTO OUTFILE "/tmp/aarAreaData.txt";

# HEADER ROW FOR SQL OUTPUT FILE
# taxon	area_no	specimen_no	ual_num	material	location	treatment	weight	flag	gain	L_Asp	D_Asp	L_Glu	D_Glu	L_Ser	D_Ser	Gly	L_Ala	hArgL	D_Ala	L_Val	D_Val	L_Phe	L_Ile	D_Phe	L_Leu	D_Alle	D_Leu	hArgL_source	c14date_no	specimen_no	labCode	ansto_num	d13c	pMC	e1Sigma	d14C	ed14C	cAgeBP	cAgeSigma	2sYng	2sOld	2sMedian	ageMedian	ageYng	ageOld	resAge	calCurve	date_prep	comments	measure2_no	user_no	specimen_no	date_measured	measured_object	x_mm	y_mm	z_mm	thick_mm	sinus_x_mm	sinus_y_mm	mass_mg	density	divide_mass_by	bias	dom	layer_no	sample_no	top	depth	bag_num	collector_no	collection_date	mass_gm	collectors	comments






# PROGRESS REPORT TO SCREEN
###################################################
print('start dataLoad-peakAreas.R')
print(paste('loading',file))

# GET - RAW DATA FILE
###################################################
areaRaw <- read.table(file, header=TRUE, na.strings=c("\\N","NA"),sep='\t')

# GET - KEY PARAMETERS / CONSTANTS
###################################################
TAXA <- levels(as.factor(areaRaw[,'taxon_no']))
TOPS <- levels(as.factor(areaRaw[,'top']))
SITES <- levels(as.factor(areaRaw[,'sample_no']))

# SET - KEY PARAMETERS / CONSTANTS
###################################################
iniA <- c('Asp','Glu','Ser','Ala','Val','Phe','Leu','A/I')
namA <- c('Aspartic','Glutamic','Serine','Alanine','Valine','Pheylalanine','Leucine','Alle/Ile')
iniA2 <- c('Asp','Glu','Ser','Ala','Val','Phe','Leu','AI')

colD <- paste('cD', iniA2,sep='_')
colL <- paste('cL', iniA2,sep='_')
colR <- paste(iniA2,'DL', sep='_')

#	CALCULATE - D & L CONCENTRATIONS
###################################################
cD_Asp <- (areaRaw[,'D_Asp']/areaRaw[,'L_hArg'])*200
cL_Asp <- (areaRaw[,'L_Asp']/areaRaw[,'L_hArg'])*200
cD_Glu <- (areaRaw[,'D_Glu']/areaRaw[,'L_hArg'])*200
cL_Glu <- (areaRaw[,'L_Glu']/areaRaw[,'L_hArg'])*200
cD_Ser <- (areaRaw[,'D_Ser']/areaRaw[,'L_hArg'])*200
cL_Ser <- (areaRaw[,'L_Ser']/areaRaw[,'L_hArg'])*200
cD_Ala <- (areaRaw[,'D_Ala']/areaRaw[,'L_hArg'])*200
cL_Ala <- (areaRaw[,'L_Ala']/areaRaw[,'L_hArg'])*200
cD_Val <- (areaRaw[,'D_Val']/areaRaw[,'L_hArg'])*200
cL_Val <- (areaRaw[,'L_Val']/areaRaw[,'L_hArg'])*200
cD_Phe <- (areaRaw[,'D_Phe']/areaRaw[,'L_hArg'])*200
cL_Phe <- (areaRaw[,'L_Phe']/areaRaw[,'L_hArg'])*200
cD_Leu <- (areaRaw[,'D_Leu']/areaRaw[,'L_hArg'])*200
cL_Leu <- (areaRaw[,'L_Leu']/areaRaw[,'L_hArg'])*200
cD_AI  <- (areaRaw[,'D_AIle']/areaRaw[,'L_hArg'])*200
cL_AI  <- (areaRaw[,'L_Ile']/areaRaw[,'L_hArg'])*200

#	CALCULATE - TOTAL CONCENTRATIONS
###################################################
cT_Asp <- areaRaw[,'D_Asp']+areaRaw[,'L_Asp']
cT_Glu <- areaRaw[,'D_Glu']+areaRaw[,'L_Glu']
cT_Ser <- areaRaw[,'D_Ser']+areaRaw[,'L_Ser']
cT_Ala <- areaRaw[,'D_Ala']+areaRaw[,'L_Ala']
cT_Val <- areaRaw[,'D_Val']+areaRaw[,'L_Val']
cT_Phe <- areaRaw[,'D_Phe']+areaRaw[,'L_Phe']
cT_Leu <- areaRaw[,'D_Leu']+areaRaw[,'L_Leu']
cT_AI <- areaRaw[,'D_AIle']+areaRaw[,'L_Ile']


#	CALCULATE - D/L RATIOS
###################################################
Asp_DL <- areaRaw[,'D_Asp']/areaRaw[,'L_Asp']
Glu_DL <- areaRaw[,'D_Glu']/areaRaw[,'L_Glu']
Ser_DL <- areaRaw[,'D_Ser']/areaRaw[,'L_Ser']
Ala_DL <- areaRaw[,'D_Ala']/areaRaw[,'L_Ala']
Val_DL <- areaRaw[,'D_Val']/areaRaw[,'L_Val']
Phe_DL <- areaRaw[,'D_Phe']/areaRaw[,'L_Phe']
Leu_DL <- areaRaw[,'D_Leu']/areaRaw[,'L_Leu']
AI_DL <- areaRaw[,'D_AIle']/areaRaw[,'L_Ile']

areaRaw <- cbind(areaRaw,cD_Asp, cL_Asp, cD_Glu, cL_Glu, cD_Ser, cL_Ser, cD_Ala, cL_Ala, cD_Val, cL_Val, cD_Phe, cL_Phe, cD_Leu, cL_Leu, cD_AI, cL_AI, cT_Asp, cT_Glu, cT_Ser, cT_Ala, cT_Val, cT_Phe, cT_Leu, cT_AI)
areaRaw <- cbind(areaRaw, Asp_DL, Glu_DL, Ser_DL, Ala_DL, Val_DL, Phe_DL, Leu_DL, AI_DL)

# DATA - SAVE CLEANED DATA
###################################################
save(areaRaw, colR, colL, colD, iniA, namA, TAXA, SITES, TOPS, file=paste(PATH.DAT,"aarAreaRaw.Rdata",sep='/'))

# PROGRESS REPORT TO SCREEN
###################################################
print('finish dataLoad-peakAreas.R')

