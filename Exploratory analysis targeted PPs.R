#Exploratory analysis of targeted polyphenol data
rawdata <- read.csv("Urinary polyphenols_EPIC.csv", header=T)
df <- rawdata[, 13:ncol(rawdata)]
pplabs <- read.csv("ppnames heatmap.csv", header=F)

#basic plots
#define and load a palette
tropical <- c("darkorange", "dodgerblue", "hotpink", "limegreen", "yellow")
palette(tropical)
#filled circles
par(pch=19)
#two columns
par(mfrow=c(2,2))
#plot scatter
plot(log(df[, 1]), col=3) #etc

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

#ECDF
library(scales)
library(ggplot2)
library(grid)
source("Base breaks function.R")

urinecor <- read.csv("urinary PPs in EPIC abbrev.csv", header=T)

#melt distributions into a single column
library(reshape2)
urinecormelt <- melt(urinecor)

#plot ecdf of each of 36 polyphenols by colour
ggplot(urinecormelt, aes(x=value, color=variable)) + stat_ecdf() +
  scale_x_continuous(name = "umol excreted per 24h", trans=log_trans(), breaks=base_breaks(), labels=prettyNum) + 
  theme_grey(base_size=10) +
  theme(legend.key.height=unit(2.5, "mm")) +
  ggsave("ECDF polyphenols2.png", height=100, width=200, units="mm")

#with latticeExtra
library(latticeExtra)
ecdfplot(~value, groups=variable, data=urinecormelt, lty=c(1:4), col=c(1:9),
         xlab="ug excreted in urine per 24 h", auto.key=TRUE,
         #main = "Excretion of 38 polyphenols in urine",
         key=list(space="right", text=list(levels(urinecormelt$variable)),
                  lines=list(lty=c(1,2,3,4), col=c(1:9))),
         scales=list(x=list(log=10, equispaced.log=F)))

key.species <- list(title="Plant", space="right", text=list(levels(CO2$Plant)),
                    points=list(pch=symbols, col=colors))

#--------------------------------------------------------------------------------------------------

#PCA by country
profilespp <- read.csv("Urinary PPs names.csv", header=T)
mat <- profilespp[, 6:41]

#impute missing values from medians for each polyphenol?
#get matrix positions of NAs
napositions <- which(is.na(mat), arr.ind=TRUE)
#insert row medians in each NA position
mat[napositions] <- apply(mat, 2, function(x) median(x, na.rm=T))[napositions[ ,2]]

#log transform and PCA
logmat <- log2(mat)
pp <- prcomp(~ ., data=logmat, scale.=T)

#recode country factor for neater labels
ctryfact <- profilespp$country
levels(ctryfact) <- c("FR","IT","GR","D")

library(ggbiplot)
ggbiplot(pp, var.axes=F, groups=profilespp$country, 
         labels=ctryfact, 
         choices=c(1,2), labels.size=3) + 
  #stat_ellipse(type="norm", size=0.5) + 
  scale_colour_brewer(palette="Set1") +
  theme_bw() + coord_fixed(ratio=0.85) + geom_vline(size=0.5, linetype="dashed") + 
  geom_hline(size=0.5, linetype="dashed") +
  theme(panel.border=element_rect(colour="black", size=0.5), panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position="none") + xlab("Principal comp. 1 (26.2% explained variation)") + 
  ylab("Principal comp. 2 (9.6% explained variation)")

ggsave("Scores targeted PPs 1v2.png", height=120, width=160, units="mm")
ggsave("Scores targeted PPs 1v2.svg", height=120, width=160, units="mm")

