tea <- read.csv("Tea new four replicates.csv", header=T, row.names=1)
#teanoimpute <- read.csv("Tea four reps no impute.csv", header=T, row.names=1)
tealabs <-c(rep("SVR", 4), rep("MK", 4), rep("BR", 4), rep("PKML", 4), rep("CF", 4),
            rep("OHO", 4), rep("DC", 4), rep("YW", 4))
logtea <- log2(tea)
#logtea2 <- log2(teanoimpute)

pclog <- prcomp(logtea, scale.=T)
#pc <- prcomp(tea, scale.=T)
pcnoscale <- prcomp(tea, scale.=F)
#pcnoimpute <- prcomp(logtea2, scale.=T)
#pcnoimpute2 <- prcomp(teanoimpute, scale.=T)

logtea <- log2(tea)
#biplot and scree plot
library(ggbiplot)
plot(pclog)
ggbiplot(pclog, var.axes=F, labels=rownames(tea), groups=tealabs) + theme_bw() +
  stat_ellipse(type="norm", linetype="dashed") +
  geom_hline(linetype="dashed") + geom_vline(linetype="dashed") +
  coord_fixed(ratio=0.75) +
  theme(panel.border=element_rect(colour="black", size=1), panel.grid.major=element_blank(),
        legend.position="none") + ggtitle("log2 transformed and UV scaled")
ggsave("Tea scores 4 reps log2 scaled2.png", height=150, width=200, units="mm")

#correlation heatmap of variables
library(corrplot)
tcor <- cor(tea, method="pearson")
corrplot(tcor, is.corr=TRUE, method="color", mar=c(1,1,1,1),
         shade.col=NA, tl.col="black", tl.srt=90, order="hclust", type="full",
         addrect=NULL, tl.cex=0.8, cl.cex=0.7)

#plot loadings separately from prcomp output
#make a data frame from the rotations calculated by prcomp
tload <- data.frame(pclog$rotation, .names = row.names(pclog$rotation))

#scatterplot to give loadings
library(ggplot2)
ggplot(tload, aes(x=PC1, y=PC2)) + geom_point() + theme_bw() + 
  geom_text(data=tload, mapping=aes(x = PC1, y = PC2+0.008, label = .names), size=4, alpha=0.8) +
  labs(x = "PC1", y = "PC2") + #ylim(-0.12, 0.08) + xlim(-0.10, 0.10) +
  geom_hline(linetype="dashed") + geom_vline(linetype="dashed") +
  theme(panel.border=element_rect(colour="black")) + ggtitle("loadings log2 transformed and UV scaled")
ggsave("Tea loadings 4 reps log2 and scaled.png", height=150, width=200, units="mm")

#cluster dendrogram
teaT <- t(tea)
hc <- hclust(dist(teaT))
plot(hc)

#ggplot2
library(reshape2)
mtcor <- data.frame(melt(tcor))
ggplot(mtcor, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + theme_bw()
