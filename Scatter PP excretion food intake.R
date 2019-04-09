#read data in and convert to matrix

food <- read.csv("Correlation PPs with foods.csv", row.names=1)
food1 <- as.matrix(food)

#melt data with reshape 
library(reshape2)
melted_food <- melt(food1)

#add two extra columns based on condition to create label groups
library(dplyr)
meltfoodcorr <- mutate(melted_food, category = ifelse(value>0.45, paste(Var2, "corr", sep="_"),
                "uncorr"), lab = ifelse(value>0.45, paste(Var1), ""))

#scatter
library(ggplot2)
library(RColorBrewer)
ggplot(meltfoodcorr, aes(x=reorder(Var2, value, FUN=median), y=value, colour=category)) +
  geom_jitter(shape=18, position=position_jitter(w=0.15, h=0.0)) +
  scale_colour_manual(values=c("orange", "brown", "darkred", "grey", "red")) +
  guides(colour=FALSE) + geom_text(aes(y=value+0.01, label=lab, hjust=1), size=2, colour="black") +
  scale_x_discrete(name="Food group") + theme_bw(base_size=10) +
  scale_y_continuous(name="Spearman rank order correlation with urinary PP", limits=c(0, 0.72)) + 
  theme(axis.text.x=element_text(angle=40, hjust=1), panel.grid.major=element_blank(),
  panel.grid.minor=element_blank()) #+ geom_hline(x=0, size=0.2)

ggsave("PP scatter correlation.png", height=100, width=130, unit="mm")

#with lattice
library(latticeExtra)
PPdot <- stripplot(value~reorder(Var1, value), data=meltfoodcorr, groups=category,
  scales=list(x=list(rot=90), y=list(tick.number=5)), ylab="Spearman correlation coefficient",
xlab="Urinary polyphenol", ylim=c(-0.02,0.75), col=c(1:6), pch=c(16:20))

#write to file
trellis.device(device="png", filename ="Dotplot PPs.png")
print(PPdot)
dev.off()

#---------------------------------------------------------------------------------------------------
#Figure for slide
#read data in and convert to matrix
library(readr)
library(tidyr)
food %>% select(-ends_with("rm")) %>% 
  gather(foodgroup, correlation, -Polyphenol) %>% ggplot(aes(x=foodgroup, y=correlation, colour=foodgroup)) + 
  scale_colour_discrete() + geom_point() + scale_x_discrete(name="Food group") + 
  scale_y_continuous(name="Spearman correlation with urinary excretion", limits=c(0, 0.72)) + 
  theme(axis.text.x=element_text(angle=40, hjust=1), legend.position="none")

ggsave("PP scatter correlation new.png", height=120, width=140, units="mm")

#--------------------------------------------------------------------------------------------------

#Scatter vs molecular weight
ue <- read.csv("Urinary polyphenols_EPIC.csv")
labs <- read.csv("uPPc codes.csv", header=T)

uelong <- ue %>% select(Center, uPPc_1:uPPc_39) %>% gather(polyphenol, conc, -Center) %>%
  filter(!is.na(conc)) %>%
  inner_join(labs, by="polyphenol") %>% 
  group_by(ppname, origin, class) %>% summarise(med=median(conc), iqr=IQR(log10(conc)))

ggplot(uelong, aes(x=iqr, y=med, colour=origin, shape=origin)) + geom_point() + #theme_bw() +
  scale_y_continuous(trans=log_trans(), breaks=base_breaks(), labels=prettyNum,
                     name=(bquote("Median urinary excretion ("*mu*"g/24h)"))) +
  #scale_shape_manual(values=c(15,16,17,25,18)) +
  scale_colour_brewer(palette="Set1") +
  #scale_colour_manual(values=c(1,2,5,4,3)) +
  scale_x_continuous(name="Variability in excretion among EPIC subjects, log10(IQR)") +
  #theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
  #legend.position=c(0.85,0.8), legend.title=element_blank())
  
  ggsave("Scatter EU variability.emf", width=160, height=120, units="mm")

library(latticeExtra)
plot(uelong$iqr, log10(uelong$med))
