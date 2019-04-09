source("Base breaks function.R")
library(tidyverse)
library(scales)

pps   <- read_csv("plasma_blk_loq3.csv") %>% gather(Compound, intensity, -ID, na.rm = T)
concs  <- pps %>% filter(ID != "Blank", ID != "LOQ")
Bl_loq <- pps %>% filter(ID == "Blank" & Compound != "Isorhamnetin" | ID == "LOQ" & Compound != "Isorhamnetin")
#Note: compounds must match in both subsets for reorder to work (removed isorhamnetin)
  
ggplot(concs, aes(reorder(Compound, intensity, FUN=median), intensity)) + 
  geom_boxplot(size=0.4, fill="grey", outlier.size = 0.5) +
  geom_point(data=Bl_loq, aes(x=Compound, y=intensity, shape=ID, colour=ID)) +
  scale_shape_manual(values=c(2,19)) + scale_colour_brewer(palette="Set1") +
  #geom_jitter(position=position_jitter(w=0.4, h=0), size=0.3, colour="grey60") + 
  coord_flip() + theme_bw() +
  scale_y_continuous(name="Plasma concentration (nM)", trans=log_trans(), 
                       breaks=base_breaks(), labels=prettyNum) +
  theme(axis.title.y=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.title=element_blank(), legend.position=c(0.9, 0.1))
  
ggsave("Boxplot plasma all subjects LOQ.svg", width=170, height=150, dpi=600, units="mm")
ggsave("Boxplot plasma all subjects LOQ2.svg")
ggsave("Boxplot plasma ICPH17 poster.svg")

tib  <- pps[1, ]

#
pps   <- read_csv("plasma_blk_loq3.csv") %>% gather(Compound, intensity, -ID, na.rm = T)
concs <- pps %>% filter(ID != "Blank", ID != "LOQ")
loq   <- pps %>% filter(ID == "LOQ" & Compound != "Isorhamnetin")
#Note: compounds must match in both subsets for reorder to work (removed isorhamnetin)

ggplot(concs, aes(reorder(Compound, intensity, FUN=median), intensity)) + 
  geom_boxplot(size=0.4, fill="grey", outlier.size = 0.5) +
  geom_point(data=loq, aes(x=Compound, y=intensity, shape=ID, colour=ID)) +
  scale_shape_manual(values=c(2,19)) + scale_colour_brewer(palette="Set1") +
  #geom_jitter(position=position_jitter(w=0.4, h=0), size=0.3, colour="grey60") + 
  coord_flip() + theme_bw() +
  scale_y_continuous(name="Plasma concentration (nM)", trans=log_trans(), 
                     breaks=base_breaks(), labels=prettyNum) +
  theme(axis.title.y=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.title=element_blank(), legend.position=c(0.9, 0.1))

