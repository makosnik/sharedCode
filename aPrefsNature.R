## NATURE TWO-COLUMN SIZE: 183mm x 247mm = 5.82 x 8.46 inches
## NATURE ONE-COLUMN SIZE: 089mm x 247mm = 2.83 x 8.46 inches

pageWidthOne <- 3.50
pageWidthTwo <- 7.20
pageHeight <- 9.72

# Helvetica, panel labels are 8pt bold (a,b,c), maximum text size 7pt, minimum text size 5pt
fontFamily <- 'Helvetica'

# lines between 0.25 and 1.0 pt

# VARIABLES ASSUMED TO BE SET BY INDIVIDUAL FILES
###################################################
#ps.options(paper="A4", horizontal=FALSE, onefile=TRUE, family="Helvetica", pointsize=14, encoding="ISOLatin1.enc")
#par(mar=c(4,4,2,2), las=0, mgp=c(2,1,0))
ps.options(paper="special", horizontal=FALSE, onefile=TRUE, family= fontFamily, pointsize=7, encoding="ISOLatin1.enc", height=pageHeight, width=pageWidthTwo)
pdf.options(paper="special", onefile=TRUE, family=fontFamily, pointsize=7, encoding="ISOLatin1.enc", height=pageHeight, width=pageWidthTwo)
par(mar=c(3,3,2,0), las=0, mgp=c(2,1,0))
