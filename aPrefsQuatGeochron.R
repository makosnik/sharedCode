## QUATERNARY GEOCHRONOLOGY TWO-COLUMN SIZE: 184mm x 240mm = 7.2440944881968 x 9.448818897648 inches
## QUATERNARY GEOCHRONOLOGY ONE-COLUMN SIZE: 089mm x 240mm = 3.5039370078778 x 9.448818897648 inches

pageWidthOne <- 3.5039370078778
pageWidthTwo <- 7.2440944881968
pageHeight <- 9.448818897648

pagePaper <- 'A4'

fontFamily <- 'Times'

# VARIABLES ASSUMED TO BE SET BY INDIVIDUAL FILES
###################################################
#ps.options(paper="A4", horizontal=FALSE, onefile=TRUE, family="Helvetica", pointsize=14, encoding="ISOLatin1.enc")
#par(mar=c(4,4,2,2), las=0, mgp=c(2,1,0))
ps.options(paper="A4", horizontal=FALSE, onefile=TRUE, family="Times", pointsize=12, encoding="ISOLatin1.enc", height=pageHeight, width=pageWidthTwo)
par(mar=c(3,3,2,0), las=0, mgp=c(2,1,0))
par(mar=c(2,4,2,1), mgp=c(2.5,0.75,0), las=1, cex=0.75, cex.lab=1.25)
