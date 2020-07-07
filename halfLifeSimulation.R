accumTotal <- 0

halfLives <- matrix(nrow=6,ncol=3)
halfLives[1,] <- c(791,935,1079)
halfLives[2,] <- c(1143,1307,1484)
halfLives[3,] <- c(453,575,704)
halfLives[4,] <- c(345,398,453)
halfLives[5,] <- c(548,667,804)
halfLives[6,] <- c(163,206,255)

rownames(halfLives) <-c('Ethalia','Marmarstoma','Pinguitellina','Liloa','Notocochlis','
Scissulina')

accumFrac <- matrix(nrow=6,ncol=3)


for (j in seq(1,6,by=1)){
for (k in seq(1,3,by=1)){
	accumTotal <- 0
	for (i in 1: accumAge){
		accumTotal <- accumTotal + accumUnit
		accumTotal <- accumTotal*(0.5^((1/halfLives[j,k])))
	}	
	accumFrac[j,k] <- accumTotal/(accumAge*accumUnit)
}
}

round(accumFrac,3)


accumTotal<- 0
accumAge <- 1000
accumUnit <- 100
halfLife <- 300

for (i in 1: accumAge){
	accumTotal <- accumTotal + accumUnit
	accumTotal <- accumTotal*(0.5^((1/halfLife)))
}

accumFrac <- accumTotal/(accumAge*accumUnit)

round((accumFrac*100),0)


(0.5^((1/halfLife)))