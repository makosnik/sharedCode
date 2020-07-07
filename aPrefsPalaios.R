## PALAIOS TWO-COLUMN SIZE: 177mm x 240mm = 7.0 x 8.5 inches
## PALAIOS ONE-COLUMN SIZE: 086mm x 240mm = 3.4 x 8.5 inches

pageWidthOne <- 3.4
pageWidthTwo <- 7.0
pageHeight <- 8.5

fontFamily <- 'Times'

# VARIABLES ASSUMED TO BE SET BY INDIVIDUAL FILES
###################################################
#ps.options(paper="A4", horizontal=FALSE, onefile=TRUE, family="Helvetica", pointsize=14, encoding="ISOLatin1.enc")
#par(mar=c(4,4,2,2), las=0, mgp=c(2,1,0))
ps.options(paper="A4", horizontal=FALSE, onefile=TRUE, family="Times", pointsize=12, encoding="ISOLatin1.enc", height=pageHeight, width=pageWidthTwo)
par(mar=c(3,3,2,0), las=0, mgp=c(2,1,0))
