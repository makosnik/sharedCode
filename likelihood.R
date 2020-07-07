#
# Burnham & Anderson 1998 p 17
# n = number of oberservations
# ml variance (var) = sum of squared deviations / n
# log likelihood = -0.5 * (n) * ln(var)
###############################################################################################
LS2LogL<-function(lmData) {

	n <- (length(summary.lm(lmData)$residuals))
	var <- (sum((summary.lm(lmData)$residuals)^2) / n )
	logl <- -0.5 * n * abs(log(var))
	
	return(logl)
}

#
# Burnham & Anderson 1998 p 48 & 51
# n = number of oberservations
# K = number of parameters
# ml variance (var) = sum of squared deviations / n
# AIC = (number of obersvations) * ln(var) + (2 * (K))
# AICc = AIC + (2K(K+1))/(n - K - 1)
#
# normal linear regression has K = 3 (intercept, slope, variance)
#
# rules of thumb (B&A p 48)
#	< 2 have substantial support
#	< 4 - 7 have considerable support
#	< 10 have essentially no support and might be omitted from further consideration.
###############################################################################################
LS2AICc<-function(lmData,K) {

	n <- (length(summary.lm(lmData)$residuals))
	var <- abs(sum(summary.lm(lmData)$residuals^2) / n )
	aicc <- (n * abs(log(var))) + (2 * K * (n / (n - K - 1)))

	return(aicc)
}