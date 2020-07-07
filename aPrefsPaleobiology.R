## PALEOBIOLOGY TWO-COLUMN SIZE: 148mm x 215mm = 5.82 x 8.46 inches
## PALEOBIOLOGY ONE-COLUMN SIZE: 072mm x 215mm = 2.83 x 8.46 inches

pageWidthOne <- 2.83
pageWidthTwo <- 5.82
pageWidthTwo <- 8.81
pageHeight <- 8.46

fontFamily <- 'Times'

# VARIABLES ASSUMED TO BE SET BY INDIVIDUAL FILES
###################################################
#ps.options(paper="A4", horizontal=FALSE, onefile=TRUE, family="Helvetica", pointsize=14, encoding="ISOLatin1.enc")
#par(mar=c(4,4,2,2), las=0, mgp=c(2,1,0))
ps.options(paper="special", horizontal=FALSE, onefile=TRUE, family="Times", pointsize=12, encoding="ISOLatin1.enc", height=pageHeight, width=pageWidthTwo)
pdf.options(paper="special", onefile=TRUE, family="Times", pointsize=12, encoding="ISOLatin1.enc", height=pageHeight, width=pageWidthTwo)
par(mar=c(3,3,2,0), las=0, mgp=c(2,1,0))
