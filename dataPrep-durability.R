###################################################
## NEEDS: dataToUse
## RETURNS: dataToUse
###################################################

print('start dataPrep-durability.R, modifies gobal dataToUse')

## TWEAK DATA
###################################################

## DEFINE DENSITY OF ARAGONITE === 2.83, CALCITE === 2.71 g/cm3
densityShell <- 2.83


###################################################
## CALCULATE SIZE [GEOMETRIC MEAN OF X*Y*Z]
###################################################

## NA - FOR DATA BELOW SCALE RESOLUTION
#for (dataCol in c('x_mm','y_mm','z_mm','thick_mm','mass_mg'))
for (dataCol in c('x_mm','y_mm','z_mm','mass_mg'))
	dataToUse[((dataToUse[,dataCol]<0.0001) & !is.na(dataToUse[,dataCol])),dataCol] <- NA

##	CALCULATE GEOMETRIC MEAN OF 3 DIMENSIONS
mmSize <- (dataToUse[,'x_mm'] * dataToUse[,'y_mm'] * dataToUse[,'z_mm'])^(1/3)
dataToUse <- cbind(dataToUse, mmSize)


###################################################
## CALCULATE VOLUME IN MM^3 [4/3*pi*r^3]
###################################################

mmVolu <- (4/3) * pi * ( (dataToUse[,'x_mm']/2) * (dataToUse[,'y_mm']/2) * (dataToUse[,'z_mm']/2) )
dataToUse <- cbind(dataToUse, mmVolu)


###################################################
## CALCULATE MASS IN MG^(1/3) [CUBE ROOT OF MASS]
###################################################

## CALCULATE MASS IN MG^(1/3)
mgMass <- (dataToUse[,'mass_mg'])
rootMass <- (mgMass)^(1/3)
dataToUse <- cbind(dataToUse,mgMass,rootMass)

## NA - FOR DATA BELOW SCALE RESOLUTION
dataToUse[((dataToUse[,'rootMass']<0.0001) & !is.na(dataToUse[,'rootMass'])),'rootMass'] <- NA


###################################################
## CALCULATE SHELL DENSITY [MASS/VOLUME (mg/mm)]
###################################################
rawDensity <- mgMass/mmVolu
relDensity <- (rawDensity / densityShell)
dataToUse <- cbind(dataToUse,rawDensity,relDensity)


###################################################
## CALCULATE DEVIATION FROM SPHERICAL
###################################################

## CALCULATE CROSS SECTIONAL AREAS [can ignore pi due to ratio]
areaYZ <- (dataToUse[,'y_mm']*dataToUse[,'z_mm'])
areaXY <- (dataToUse[,'x_mm']*dataToUse[,'y_mm'])
areaXZ <- (dataToUse[,'x_mm']*dataToUse[,'z_mm'])
areaSP <- (dataToUse[,'mmSize']*dataToUse[,'mmSize'])
dataToUse <- cbind(dataToUse,areaYZ,areaXY,areaXZ,areaSP)

areaMin <- vector(length=nrow(dataToUse))
## WHAT I WANT? IS THE SMALLEST CROSS SECTIONAL AREA AS THIS SHOULD BE MOST LIKELY TO BREAK?
for (i in 1:nrow(dataToUse)) {
	if ((!is.na(areaYZ[i])) && (!is.na(areaXY[i])) && (!is.na(areaXZ[i]))) {
		
		if ((areaYZ[i] < areaXY[i]) && (areaYZ[i] < areaXZ[i])) {
			areaMin[i] <- areaYZ[i]
		} else {
			if ((areaXY[i] < areaYZ[i]) && (areaXY[i] < areaXZ[i])) {
				areaMin[i] <- areaXY[i]
			} else {
				areaMin[i] <- areaXZ[i]
			}
		}
	} else {
		areaMin[i] <- NA
	}
}
dataToUse <- cbind(dataToUse,areaMin)

relXsection <- areaMin/areaSP
dataToUse <- cbind(dataToUse, relXsection)


###################################################
##	OTHER "CONSTANTS"
###################################################

## DEFINE MASS OF SPHERE OF SET SIZE
sizeMax <- ceiling(max(dataToUse[,'mmSize'],na.rm=TRUE))
sizRange <- c(0.1,sizeMax)
volRange <- (4/3)*pi*((sizRange/2)^3)
weiRange <- (volRange*densityShell)^(1/3)
w <- length(weiRange)

##	MINIMUM / MAXIMUM MASS
massMax <- ceiling(max((dataToUse[,'rootMass']),na.rm=TRUE))
massMin <- floor(min((dataToUse[,'rootMass']),na.rm=TRUE))

##	MINIMUM / MAXIMUM DENSITY
densMax <- ceiling(max((dataToUse[,'relDensity']),na.rm=TRUE))
densMin <- floor(min((dataToUse[,'relDensity']),na.rm=TRUE))

##	MINIMUM / MAXIMUM THICKNESS
#thickMax <- ceiling(max(dataToUse[,'thick_mm']/dataToUse[,'mmSize'],na.rm=TRUE)*10)/10
#thickMin <- floor(min(dataToUse[,'thick_mm']/dataToUse[,'mmSize'],na.rm=TRUE)*10)/10

##	MINIMUM / MAXIMUM DEVIATION
devMax <- ceiling(max((dataToUse[,'areaMin']/dataToUse[,'areaSP']),na.rm=TRUE)*10)/10
devMin <- floor(min((dataToUse[,'areaMin']/dataToUse[,'areaSP']),na.rm=TRUE)*10)/10

xla <- expression(paste('size ',sqrt('xyz',3),' (mm)'))

