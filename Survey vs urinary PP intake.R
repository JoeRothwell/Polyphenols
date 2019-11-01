hdr <- read.csv("data/HDR.csv", row.names=1)
ffq <- read.csv("data/FFQ.csv", row.names=1)
hdr1 <- read.csv("data/HDRaglyc.csv", row.names=1)
ffq1 <- read.csv("data/FFQaglyc.csv", row.names=1)

#read in urinary polyphenol names in data frame order
ppnames <- scan("data/ppnames heatmap.csv", what="character", sep=",")

library(gplots)
#for 24HDR and FFQ
heatmap.2(as.matrix(hdr), col=redblue(256), trace="none", labCol = ppnames, 
          margins=c(10,2), offsetCol = 0.02, offsetRow = 0.1, cexRow = 0.1, cexCol = 0.6)

heatmap.2(as.matrix(ffq), col=redblue(256), trace="none", labCol = ppnames,
          margins=c(10,2), offsetCol = 0.02, offsetRow = 0.1, cexRow = 0.1, cexCol = 0.6)

ppnames <- scan("ppnames heatmap.csv", what="character", sep=",")

#for 24HDR and FFQ aglycones
heatmap.2(as.matrix(hdr1), col=redblue(256), trace="none", labCol = ppnames,
          margins=c(10,4), offsetCol = 0.02, offsetRow = 0.1, cexRow = 0.2, cexCol = 0.6)

heatmap.2(as.matrix(ffq1), col=redblue(256), trace="none", labCol = ppnames,
          margins=c(10,5), offsetCol = 0.02, offsetRow = 0.1, cexRow = 0.2, cexCol = 0.6)