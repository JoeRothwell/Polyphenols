library(ggplot2)
library(dplyr)

tea <- read.csv("Tea new four replicates.csv", header=T, row.names=1)
tealabs <-c(rep("SVR", 4), rep("MK", 4), rep("BR", 4), rep("PKML", 4), rep("CF", 4),
            rep("OHO", 4), rep("DC", 4), rep("YW", 4))
logtea <- log2(tea)

pclog <- prcomp(logtea, scale.=T)
pcnoscale <- prcomp(tea, scale.=F)

groups <- c(rep("SVR", 4), rep("MK", 4), rep("BR", 4), rep("PKML", 4), rep("CF", 4),
            rep("OHO", 4), rep("DC", 4), rep("YW", 4))

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
