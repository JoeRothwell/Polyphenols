tea <- read.csv("Tea_PCA_pp.csv", header=T)
#new set of data 19-feb 2015
teanew <- read.csv("Teanew.csv", header=T)
#discard altitude data and subset observation names separately
cmpds <- teanew[ , 5:71]

pc <- prcomp(cmpds, .scale=FALSE)
tlabs <- teanew[, 1]
groups <- teanew[ , 4]
#biplot and scree plot
biplot(pc, xlabs=tlabs)
plot(pc)

#correlation heatmap of variables
library(corrplot)
tcor <- cor(cmpds, method="spearman")
corrplot(tcor, is.corr=TRUE, method="color", mar=c(1,1,1,1),
         shade.col=NA, tl.col="black", tl.srt=90, order="hclust", type="full",
         addrect=NULL, tl.cex=0.8, cl.cex=0.7)

library(ggbiplot)
ggbiplot(pc, choices = c(1,2), labels=tlabs, var.axes=F) + theme_bw() + stat_ellipse(type="norm", size=0.2) +
  geom_hline(yintercept=0, linetype="solid", size=0.2) + geom_vline(xintercept=0, linetype="solid", size=0.2) +
  theme(panel.border=element_rect(colour="black", size=1), panel.grid.major=element_blank())
ggsave("Tea scores.png", height=150, width=150, units="mm") 

#plot loadings separately from prcomp output
#make a data frame from the rotations calculated by prcomp
tload <- data.frame(pc$rotation, .names = row.names(pc$rotation))

#scatterplot to give loadings
library(ggplot2)
ggplot(tload, aes(x=PC1, y=PC2)) + geom_point() + theme_bw() + 
  geom_text(data=tload, mapping=aes(x = PC1, y = PC2+0.015, label = .names), size=4, alpha=0.8) +
  labs(x = "PC1", y = "PC2") + #ylim(-0.12, 0.08) + xlim(-0.10, 0.10) +
  geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed") +
  theme(panel.border=element_rect(colour="black"))
ggsave("Tea loadings2.png", height=150, width=200, units="mm")

#cluster dendrogram
cmpdsT <- t(cmpds)
hc <- hclust(dist(cmpdsT))
plot(hc)

#ggplot2
library(reshape2)
mtcor <- data.frame(melt(tcor))
ggplot(mtcor, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + theme_bw()

