#Exploratory analysis of targeted polyphenol data
rawdata <- read.csv("data/Urinary polyphenols_EPIC.csv", header=T)
df <- rawdata[, 13:ncol(rawdata)]
pplabs <- read.csv("data/ppnames heatmap.csv", header=F)


#basic plots
#define and load a palette
tropical <- c("darkorange", "dodgerblue", "hotpink", "limegreen", "yellow")
palette(tropical)
#filled circles
par(pch=19)
#two columns
par(mfrow=c(2,2))
#plot scatter
plot(log(df[, 1]), col=3)
plot(log(df[, 2]), col=2)
plot(log(df[, 3]), col=1)
plot(log(df[, 4]), col=4)

#boxplot
boxplot(log(df), horizontal=T, col=5)
boxplot(log(df[, 4]), col=2)

#histogram
par(mfrow=c(1,2))
hist(log(df[, 3]), col=2, breaks=50)
hist(log(df[, 4]), col=3, breaks=50)

plot(density(df[, 4]), col=2)
#to overlay another
lines(density(df[, 9]), col=3)

library(gplots)
colramp = colorRampPalette(c(2, "white", 3))(15)
heatmap.2(log(as.matrix(df)), scale="column", trace="none", labCol = unlist(pplabs[-c(9,10)]),
          offsetCol=0.1)

# PCPR2
# Subset compound concentrations and metadata
library(pcpr2)
library(zoo)
library(tidyverse)

mat <- rawdata %>% select(starts_with("uPPc"))
X_DataMatrix <- na.aggregate(mat, FUN = function(x) min(x)/2) %>% scale

# Y-variables. Remove country as directly correlated with centre
Z_InterestFactors <- rawdata %>% select(Batch:RE_ENERGY, -Country)
dev.off()
plotProp(runPCPR2(X_DataMatrix, Z_InterestFactors))
