
# BOOTSTRAP FUNCTIONS
###################################################

# SAMPLE MEDIAN
###################################################
bootMedian <- function(x, d) { 
  return(median(x[d])) 
} 

# bootFit <- boot(datesTaxon[,colAges[age]], samplemedian, R=1000)

# BOOTSTRAP EXPONENTIAL DISTRIBUTION FIT
#
# 		bootFit <- boot(histData+MIN, bootFitExponential, R=REPS)
#		bootCI <- boot.ci(bootFit, type="perc")
#
###################################################
bootFitExponential <- function(x, d) {
	return(fitdistr(x[d], densfun='exponential')$estimate)
}

bootFitGeometric <- function(x, d) {
	return(fitdistr(x[d], densfun='geometric')$estimate)
}

# BOOTSTRAP OPIMAL REGRESSION
###################################################
bootRegressOpim <- function(x, y, d) {
	return(regressOpim(x[d], y[d]))
}

# BOOTSTRAP CALIBRATION CURVE
###################################################
bootCalCurveOptimExponent <- function(x, d, w) {
	return(aaCalCurveOptimExponent(c(0,6), data14C=x[d,c(1:3)], dataAmino=x[d,4], weight=w))
}