# Exploratory analysis of targeted polyphenol data
rawdata <- read.csv("data/Urinary polyphenols_EPIC.csv", header=T)
df <- rawdata[, 13:ncol(rawdata)]
pplabs <- read.csv("data/ppnames heatmap.csv", header=F)

# basic plots: define palette
tropical <- c("darkorange", "dodgerblue", "hotpink", "limegreen", "yellow")
palette(tropical)
#filled circles, 2 columns
par(pch=19)
par(mfrow=c(2,2))
#plot scatter for cols 1 to 4
for(i in 1:4) plot(log(df[, i]), col = i)

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
library(pcpr2)
library(tidyverse)
mat <- rawdata %>% select(starts_with("uPPc"))
X_DataMatrix <- zoo::na.aggregate(mat, FUN = function(x) min(x)/2) %>% scale

# Y-variables. Remove country as directly correlated with centre
Z_InterestFactors <- rawdata %>% select(Batch:RE_ENERGY, -Country)
dev.off()
output <- runPCPR2(X_DataMatrix, Z_InterestFactors)
plot(output, main = "34 polyphenols")
