source('../aarCalibration/calibrationFunctions.R')
source('../sharedCode/regress.R')


dataAll <- dataToUse

taxon <- 'Liloa'
dataAll <- dataAll[(dataAll[,'material']==taxon),]
amino <- 1


# SET AD XXXX TO 0 NOT AD 1950
dateMod <- 56
DATEZERO <- 1950 + dateMod

dataAll[,'X2sMedian'] <- dataAll[,'X2sMedian'] + dateMod
dataAll[,'X2sYng'] <- dataAll[,'X2sYng'] + dateMod
dataAll[,'X2sOld'] <- dataAll[,'X2sOld'] + dateMod

#dataCal <- dataAll[!is.na(dataAll[,'X2sMedian']),]

ageMin <- -100
ageMax <- max(datesToUse[,'X2sOld'],na.rm=TRUE)


w1 <- aaCalCurveOptim(c(0.4,6), ageMin, ageMax, data14C= dataAll[,c('X2sMedian','X2sYng','X2sOld')], dataAmino= dataAll[,colR[amino]], nameAmino=iniA[amino], nameSpecimen=taxon, DATEZERO=DATEZERO, LIVING=LIVING)

dataAll <- cbind(dataAll,w1)