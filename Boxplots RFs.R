library(ggplot2)

#Will's boxplot

g <- ggplot(boxplot.df, aes(x=Box, y=XObs_gLog_QC.LSC)) +
  geom_boxplot(outlier.colour="red", outlier.size=0) +
  stat_summary(fun.y=mean, geom="point", shape=18, size=12, col="red", fill="red") +
  geom_jitter(position=position_jitter(w=0.15,h=0.15)) +
  geom_hline(yintercept=b1$stats[[5]], colour="red", size=1) +
  annotate("text", label=Above.zero.samples, y=max(nonzero.boxplot[,j])*1.1, x=2, size=3.5) +
  annotate("text", label=Zero.outliers, y=max(nonzero.boxplot[,j])*1.1,x=1,size=3.5) +
  labs(title=Titlename) + theme_bw(20)

#Box and points 16 processes

rfs <- read.csv("raw_data.csv", header=T)

fig3 <- ggplot(rfs, aes(x=process.full, y=mean.rf)) +
  geom_boxplot(outlier.colour="red",outlier.size=0) +
  geom_jitter(position=position_jitter(w=0.15, h=0), size=2, colour="grey60", alpha=.6) +
  #stat_summary(fun.y=mean, geom="point",shape=18, size=4) +
  geom_hline(yintercept=1, linetype="dashed") + coord_flip() +
  theme_bw(base_size=9) + 
  theme(panel.border=element_rect(colour="black"), axis.title.y=element_blank(),
        panel.grid.major=element_blank()) +
  scale_y_continuous(name="Retention factor", trans=log_trans(), breaks=base_breaks(), labels=prettyNum) +
  scale_x_discrete(limits=c("Pressure-boiled", "Pressure-steamed", "Microwaved", "Canned",
                            "Steamed", "Boiled", "Fried", "Frozen", "Blanched", "Jam making",
                            "Stored at room temperature", "Frozen, stored frozen",
                            "Stored refrigerated", "Stored frozen", "Pasteurized", "Grilled"),
                   labels=c("Pressure\nboiled", "Pressure\nsteamed", "Microwaved", "Canned",
                            "Steamed", "Boiled", "Fried", "Frozen", "Blanched", "Jam-\nmaking",
                            "Stored at\nroom temp", "Frozen,\nstored frozen",
                            "Stored\nrefrigerated", "Stored\nfrozen", "Pasteurized", "Grilled"))

ggsave("fig3revised.png", width=80, height=110, units="mm", dpi=600, plot=fig3)

#in lattice
library(lattice)
bwplot(reorder(process, mean.rf, FUN=median) ~ mean.rf, data=rfs, jitter.data=T,
       xlab="Mean retention factor",
       scales=list(x=list(log=10, equispaced.log=F)), subset=rfs$top.process=="yes")

#----------------------------------------------------------------------------------------------------

#Boxplot all RFs by method
ggplot(rfs, aes(x=method, y=mean.rf)) + geom_boxplot(outlier.size=1.5, varwidth=T) + theme_bw() +
  geom_hline(yintercept=1, linetype="dashed") + theme(axis.title.x=element_blank()) +
  scale_y_log10(name="Retention factor", breaks=c(0.01,0.05,0.1,0.5,1,5,10,50)) + 
  scale_x_discrete(labels=c("Chromatography", "Chromatography\nwith\nhydrolysis",
                            "Folin-\nCiocalteu\nassay", "pH\ndifferential\nassay"))

ggsave("Boxplot all RFs by method.png", width=133, height=100, unit="mm")

#with base graphics and nice colours
library(RColorBrewer)
boxplot(log(rfs$mean.rf) ~ method, data = rfs, xlab="Retention factor", horizontal = T, col=brewer.pal(4,"Set3"))

#Frequency polygon
ggplot(rfs, aes(x=nraw.rf)) + geom_freqpoly(colour="black") + theme_bw() + 
  labs(x="Raw retention factors", y="Aggregated retention factors") 
ggsave("aggregation histogram.png", width=133, height=100, unit="mm")

#Boxplot RFs by method
ggplot(rfs, aes(x=method, y=mean.rf)) + geom_boxplot(outlier.size=1.5) + theme_bw() +
  geom_hline(yintercept=1, linetype="dashed") + labs(x="", y="Retention factor") +  
  scale_y_continuous(expand = c(0.02,0.02), limits=c(0,7.5), breaks=c(0:8))

ggsave("Boxplot all RFs by method.png", width=4, height=3, unit="in")

