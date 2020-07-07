# CALCULATE THE DISTANCE BETWEEN TWO POINTS
###################################################
distancePythagorean <- function(a,b) {

	dA <- (a[1] - b[1])
	dB <- (a[2] - b[2])
	dC <- sqrt(dA^2 + dB^2)

	return (dC)
}

# CALCULATE THE DISTANCE BETWEEN A SITES AND A LIST OF POINTS
###################################################
distanceToAllPoints <- function(point,objects) {
	
	nObjects <- length(objects[,1])
	
	distanceVector <- vector(mode = "numeric", length = nObjects)

	for (i in 1:nObjects)
		distanceVector[i] <- distancePythagorean(point,objects[i,])
		
	return(distanceVector)
}

# CALCULATE THE DISTANCES BETWEEN A LIST OF SITES AND A LIST OF POINTS
###################################################
distancesAll2All <- function(sites,objects) {

	nSites <- length(sites[,1])
	nObjects <- length(objects[,1])

	distanceMatrix <- matrix(data=0, nrow = nSites, ncol = nObjects)

	for (i in 1:nSites)
		distanceMatrix[i,] <- distanceToAllPoints(sites[i,],objects)
		
	return(distanceMatrix)

}

# CALCULATE THE DISTANCES BETWEEN A LIST OF SITES AND A LIST OF POINTS
###################################################
distancesSelected <- function(distances,number) {

	nSites <- length(distances[,1])

	distanceVector <- vector(mode = "numeric", length = nSites)

	for (i in 1:nSites) {
		dVec1 <- sort(distances[i,])
		distanceVector[i] <- median(dVec1[1:number])
	}
	
	return(distanceVector)
}