#boxplot with or without points
#prepare data by subsetting and melting individual and loq and blank data
library(reshape2)
stackpp <- melt(pps[1:475, ], id.vars="ID", variable.name="polyphenol", value.name="concentration", na.rm=T)
#stackloqblank <- melt(pps[476:477, ], id.vars="ID", variable.name="polyphenol", value.name="concentration")
stackloq <- melt(pps[476, ], id.vars="ID", variable.name="polyphenol", value.name="concentration")

#join codes and polyphenol names together
library(dplyr)
stackpplabs <- inner_join(stackpp, labs, by="polyphenol")
stackloqlabs <- inner_join(stackloq, labs, by="polyphenol")

#plot overlaying the boxplot geom with another point geom
library(ggplot2)
library(scales)
source("Base breaks function.R")
ggplot(stackpplabs, aes(x=reorder(ppname, concentration, FUN=median), y=concentration)) + theme_bw() + 
  #geom_boxplot(outlier.size=0) + 
  geom_jitter(position=position_jitter(w=0.2, h=0), size=1, colour="grey60", alpha=.6) +
  geom_boxplot(outlier.size=0, alpha=.2, colour="blue", fill="blue", size=0.4) + 
  labs(x="Polyphenol", y="Urinary excretion (um)") +
  geom_point(data=stackloqlabs, aes(x=ppname, y=concentration, shape=ID, colour=ID)) +
  coord_flip() + scale_shape_manual(values=23) +
  scale_colour_manual(values=10) +
  theme(legend.title=element_blank(),
        axis.title.y=element_blank(),
        #axis.text.x=element_text(angle=45,  hjust=1,  vjust=1),
        #panel.grid.major=element_blank(),
        #legend.position=c(0.1, 0.85)) +
        legend.position=c(0.85, 0.1)) +
  scale_y_continuous(trans=log_trans(), breaks=base_breaks(), labels=prettyNum) +
  ggsave("Boxplot PPs manuscript.png", width=140, height=130, units="mm")

#----------------------------------------------------------------------------------------------------

#boxplot without points, percentile 10:90
#prepare data by subsetting and melting individual and loq and blank data
stackpp <- melt(pps[1:475, ], id.vars="ID", variable.name="polyphenol", value.name="concentration", na.rm=T)
stackloqblank <- melt(pps[476:477, ], id.vars="ID", variable.name="polyphenol", value.name="concentration")

#function to give boxplot whiskers to the 10th and 90th percentiles
f <- function(x) {
  r <- quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

#plot overlaying the boxplot geom with another point geom
ggplot(stackpp, aes(x=reorder(polyphenol, concentration, FUN=median), y=concentration)) + theme_bw() + 
  stat_summary(fun.data = f, geom="boxplot") +
  #geom_boxplot(outlier.size=1) +
  labs(x="Polyphenol", y="Urinary concentration (um)") +
  geom_point(data=stackloqblank, aes(x=polyphenol, y=concentration, shape=ID, colour=ID)) +
  scale_shape_manual(values=c(1, 3)) + #coord_flip() +
  scale_colour_brewer(palette="Set1") +
  theme(legend.title=element_blank(),
        axis.text.x=element_text(angle=45,  hjust=1,  vjust=1),
        panel.grid.major=element_blank(),
        legend.position=c(0.1, 0.85)) +
  #legend.position=c(0.85, 0.1)) +
  scale_y_continuous(trans=log_trans(), breaks=base_breaks(), labels=prettyNum) +
  ggsave("Calibration PPs 1090.png", width=200, height=150, units="mm")

#library(latticeExtra)
#bwplot(concentration ~ polyphenol, data=stackpp, xlab="Polyphenol", horizontal=F)