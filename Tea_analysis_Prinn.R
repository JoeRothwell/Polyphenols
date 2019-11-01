# Tea analysis for Prinn

# Original dataset ----
# had one replicate per sample only
# tea <- read.csv("data/Tea_PCA_pp.csv", header=T)
# new set of data 19-feb 2015
teanew <- read.csv("data/Teanew.csv", header=T)

# discard altitude data and subset observation names separately
cmpds <- teanew[ , 5:71]
pc <- prcomp(cmpds, .scale=FALSE)
tlabs <- teanew[, 1]
groups <- teanew[ , 4]
#biplot and scree plot
biplot(pc, xlabs=tlabs)
plot(pc)

# Spearman correlation heatmap of variables
library(corrplot)
tcor <- cor(cmpds, method="spearman")
corrplot(tcor, is.corr=TRUE, method="color", mar=c(1,1,1,1),
         shade.col=NA, tl.col="black", tl.srt=90, order="hclust", type="full",
         addrect=NULL, tl.cex=0.8, cl.cex=0.7)

# Scatterplot to give loadings
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

# New data: four replicates per sample ----
tea <- read.csv("data/Tea new four replicates.csv", header=T, row.names=1)
tea0 <- read.csv("data/Tea four reps no impute.csv", header=T, row.names=1)

tealabs <- rep(c("SVR", "MK", "BR", "PKML", "CF", "OHO", "DC", "YW"), each = 4)

logtea <- log2(tea)
#logtea2 <- log2(tea0)

pclog <- prcomp(logtea, scale.=T)
#pc <- prcomp(tea, scale.=T)
pcnoscale <- prcomp(tea, scale.=F)
#pcnoimpute <- prcomp(logtea2, scale.=T)
#pcnoimpute2 <- prcomp(tea0, scale.=T)

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
  geom_hline(yintercept = 0, linetype="dashed") + geom_vline(xintercept = 0, linetype="dashed") +
  theme(panel.border=element_rect(colour="black")) + 
  ggtitle("loadings log2 transformed and UV scaled")
ggsave("Tea loadings 4 reps log2 and scaled.png", height=150, width=200, units="mm")

#cluster dendrogram
teaT <- t(tea)
hc <- hclust(dist(teaT))
plot(hc)

#ggplot2
library(reshape2)
mtcor <- data.frame(melt(tcor))
ggplot(mtcor, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + theme_bw()

# New PCA (formerly Tea new PCA.R)
pclog <- prcomp(logtea, scale.=T)
pcnoscale <- prcomp(tea, scale.=F)
groups <- tealabs

d <- data.frame(pclog$x) %>% mutate(Tea = groups)
e <- data.frame(unscaled$x) %>% mutate(Tea = groups)

#ggplot2 settings
theme_pca <- theme_bw() + theme(panel.border=element_rect(colour="black", size=1), 
                                panel.grid.major=element_blank(), legend.key = element_rect(colour="white"), legend.title=element_blank())

ggplot(d, aes(x=PC1, y=PC2)) + geom_point(data=d, mapping = aes(shape=Tea)) + 
  stat_ellipse(type="norm", size=0.2) + geom_hline(yintercept = 0, size=0.2) + theme_pca +
  geom_vline(xintercept = 0, size=0.2) + xlab("Scores on PC1") + ylab("Scores on PC2") +
  scale_shape_manual(values=c(1,2,5,6,15,17,18,19)) + xlab("Scores on PC1") + ylab("Scores on PC2")

ggplot(e, aes(x=PC1, y=PC2)) + geom_point(data=e, mapping = aes(shape=Tea)) + 
  stat_ellipse(type="norm", size=0.2) + theme_pca + geom_hline(yintercept = 0, size=0.2) + 
  geom_vline(xintercept = 0, size=0.2) +
  scale_shape_manual(values=c(1,2,5,6, 15, 17, 18,19)) + xlab("Scores on PC1") + ylab("Scores on PC2")

ggsave("Scores plot symbols scaled bw.png", height=100, width=160, units="mm")

#make a data frame from the rotations calculated by prcomp
f <- data.frame(pclog$rotation, .names = c(1:nrow(pclog$rotation)))
g <- data.frame(pcnoscale$rotation, .names = c(1:nrow(pcnoscale$rotation)))

#scatterplot to give loadings
ggplot(f, aes(x=PC1, y=PC2)) + geom_point() + theme_pca + 
  geom_hline(yintercept = 0, size=0.2) + theme_pca + geom_vline(xintercept = 0, size=0.2) +
  geom_text(data=f, mapping=aes(x = PC1-0.002, y = PC2-0.009, label = .names), size=3, alpha=0.8) +
  labs(x = "Loadings on PC1", y = "Loadings on PC2")

ggplot(g, aes(x=PC1, y=PC2)) + geom_point() + theme_pca + 
  geom_hline(yintercept = 0, size=0.2) + theme_pca + geom_vline(xintercept = 0, size=0.2) +
  geom_text(data=g, mapping=aes(x = PC1-0.005, y = PC2-0.025, label = .names), size=3, alpha=0.8) +
  labs(x = "Loadings on PC1", y = "Loadings on PC2")

ggsave("Loadings plot numbers unscaled.png", height=100, width=160, units="mm")

#filtered variables
filter <- c(3,6,11,12,14,18,19,29,31,34,45) #positions of variables to be retained (from email 22/2/16)
teafilter <- logtea[, filter]
pclog <- prcomp(teafilter, scale.=T)
unscaled <- prcomp(teafilter, scale.=F) #after this obtain d and e variables from above

#ggplot2 settings
theme_pca <- theme_bw() + theme(panel.border=element_rect(colour="black", size=1), 
                                panel.grid.major=element_blank(), legend.key = element_rect(colour="white"), legend.title=element_blank())

ggplot(e, aes(x=PC1, y=PC2)) + geom_point(data=e, mapping = aes(shape=Tea)) + 
  stat_ellipse(type="norm", size=0.2) + geom_hline(yintercept = 0, size=0.2) + theme_pca +
  geom_vline(xintercept = 0, size=0.2) + xlab("Scores on PC1") + ylab("Scores on PC2") +
  scale_shape_manual(values=c(1,2,5,6,15,17,18,19)) + xlab("Scores on PC1") + ylab("Scores on PC2")  

ggsave("Scores reduced variables unscaled.png", height=100, width=160, units="mm")

# PCPR2 (formerly PCPR2 Tea.R)

tea <- read.csv("data/Tea mod joe.csv")

# Subset X and Y variables
X_DataMatrix <- log2(tea[, 6:72]) %>% as.matrix
Y_var <- tea[, 2:4]
Y_var <- Y_var %>% mutate_at(vars(Replicate, Experiment), as.factor)

library(pcpr2)
output <- runPCPR2(X_DataMatrix, Y_var)
plot(output, col = "red")

