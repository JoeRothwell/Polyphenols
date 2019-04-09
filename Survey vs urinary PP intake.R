hdr <- read.csv("HDR.csv", header=T, row.names=1)
ffq <- read.csv("FFQ.csv", header=T, row.names=1)
hdraglyc <- read.csv("HDRaglyc.csv", header=T, row.names=1)
ffqaglyc <- read.csv("FFQaglyc.csv", header=T, row.names=1)

#read in urinary polyphenol names in data frame order
ppnames <- scan("ppnames heatmap.csv", what="character", sep=",")

library(gplots)
#for 24HDR
heatmap.2(as.matrix(hdr), col=redblue(256), trace="none", labCol = ppnames, 
          margins=c(10,2), offsetCol = 0.02, offsetRow = 0.1, cexRow = 0.1, cexCol = 0.6)

#for FFQ
heatmap.2(as.matrix(ffq), col=redblue(256), trace="none", labCol = ppnames,
          margins=c(10,2), offsetCol = 0.02, offsetRow = 0.1, cexRow = 0.1, cexCol = 0.6)

ppnames <- scan("ppnames heatmap.csv", what="character", sep=",")

#for 24HDR aglycones
heatmap.2(as.matrix(hdraglyc), col=redblue(256), trace="none", labCol = ppnames,
          margins=c(10,4), offsetCol = 0.02, offsetRow = 0.1, cexRow = 0.2, cexCol = 0.6)

#for FFQ aglycones
heatmap.2(as.matrix(ffqaglyc), col=redblue(256), trace="none", labCol = ppnames,
          margins=c(10,5), offsetCol = 0.02, offsetRow = 0.1, cexRow = 0.2, cexCol = 0.6)