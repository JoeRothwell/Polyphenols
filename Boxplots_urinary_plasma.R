source("Base breaks function.R")
library(tidyverse)
library(ggplot2)
library(scales)

ppwide <- read.csv("Urinary PPs LOQ blank uPPc.csv", header=T)
labs   <- read.csv("uPPc codes.csv", header=T)

#gather all variables except batch and ID and call the new columns polyphenol and conc. Must remove NAs
pplong <- gather(ppwide, polyphenol, conc, -Batch, -ID, na.rm=T)
pplabs <- pplong %>% inner_join(labs, by="polyphenol") %>% filter(Batch != 99)
loqs   <- pplong %>% inner_join(labs, by="polyphenol") %>% filter(ID == "LOQ")

ggplot(pplabs, aes(x=reorder(ppname, conc, FUN=median), y=conc)) + theme_bw() + coord_flip() +
  geom_boxplot(outlier.size=0, size=0.4) +
  #add the individual points:
  geom_jitter(position=position_jitter(w=0.4, h=0), size=0.2, colour="grey60") +
  #add the LOQ point for each
  geom_point(data=loqs, aes(x=ppname, y=conc, shape=ID, colour=ID)) +
  scale_shape_manual(values=1) + scale_colour_brewer(palette="Set1") +
  scale_y_continuous(trans=log_trans(), breaks=base_breaks(), labels=prettyNum,
                     name=(bquote("Urinary concentration ("*mu*"M)"))) +
  theme(legend.title=element_blank(), axis.title.y=element_blank(),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position=c(0.85, 0.1))
ggsave("Boxplot poster David points.png", width=170, height=190, dpi=600, units="mm")
ggsave("Boxplot poster David points.svg", width=170, height=190, units="mm")

#ID != "35____35814594"

#same for plasma data for 200 subjects 16 May 2016  -----------------------------------------------

plasma <- read.csv("EPIC plasma 200 subjects.csv", header=T)
pplong <- plasma %>% select(-uPPc_30, -uPPc_10) %>% filter(ID != "35____35814594") %>%
          gather(polyphenol, conc, -Batch, -ID, na.rm=T)
pplabs <- pplong %>% inner_join(labs, by="polyphenol") %>% filter(Batch != 99)
loqs   <- pplong %>% inner_join(labs, by="polyphenol") %>% filter(Batch == 99)

ggplot(pplabs, aes(x=reorder(ppname, conc, FUN=median), y=conc)) + theme_bw() + coord_flip() +
  geom_boxplot(outlier.size = 0, size=0.4) +
  geom_jitter(position=position_jitter(w=0.4, h=0), size=0.3, colour="grey60") +
  geom_point(data=loqs, aes(x=ppname, y=conc, shape=ID, colour=ID)) +
  scale_shape_manual(values=1) + scale_colour_brewer(palette="Set1") +
  scale_y_continuous(name="Plasma concentration (nM)", trans=log_trans(), 
                     breaks=base_breaks(), labels=prettyNum) +
  theme(axis.title.y=element_blank(), panel.grid.major=element_blank(),legend.title=element_blank(),
legend.position=c(0.85, 0.1), panel.grid.minor=element_blank())
ggsave("Boxplot EPIC plasma 200 subjects.png", width=170, height=190, dpi=600, units="mm")

#generate correlation matrix
rownames(plasma) <- NULL
tcor <- cor(plasma[-114, -c(1,11,30)], use="pairwise.complete.obs", method="pearson")

library(corrplot)
corrplot(tcor, mar=c(1,1,1,1), is.corr=FALSE, method="color",
         shade.col=NA, tl.col="black", tl.srt=90, order="hclust", type="full",
         tl.cex=0.8, cl.cex=0.7)

# Segplot for poster ----

library(latticeExtra)

# Example
data(USCancerRates)
segplot(reorder(factor(county), rate.male) ~ LCL95.male + UCL95.male,
        data = subset(USCancerRates, state == "Washington"),
        draw.bands = FALSE, centers = rate.male)

# Polyphenol concentrations
ue   <- read.csv("Urinary polyphenols_EPIC.csv")
labs <- read.csv("uPPc codes.csv", header=T)

uelong <- ue %>% select(Center, uPPc_1:uPPc_39) %>% gather(polyphenol, conc, -Center) %>%
  filter(!is.na(conc)) %>% inner_join(labs, by="polyphenol") %>% 
  group_by(ppname, origin, class) %>% summarise(med=median(conc), low=quantile(conc, 0.25),
                                                high=quantile(conc, 0.75))

#reorder factor levels
uelong$origin <- factor(uelong$origin, levels=c("Food", "Endogenous", "Microbiota", 
                                                "Food/Endogenous", "Food/Microbiota"))

#solution to colour problem found at 
#http://stackoverflow.com/questions/27435642/different-interval-colors-in-a-segplot-in-r
segpp <- segplot(reorder(ppname, med) ~ log10(low) + log10(high), data=uelong, draw.bands=F,
                centers = log10(med), xlab=(bquote("Median excretion ("*mu*"g/24h)")), groups=origin,
                scales = list(x = list(labels=c(0.01, 0.1, 1, 10, 100)), tick.number=5),
                panel = function(...) { #panel.abline(v=1, lty=3)
                  panel.superpose(...)
                },
                panel.groups = function (x,y,z,subscripts,col,col.line,centers,...){
                panel.segplot(x,y,z[subscripts], centers=centers[subscripts], subscripts=T,
                  col=col.line,...)
                })

#write to file
trellis.device(device="svg", filename ="Segplot PPs.svg")
print(segpp)
dev.off()

# All plasma concentrations (formerly Boxplot David all plasma concs.R)
pps   <- read_csv("data/plasma_blk_loq3.csv") %>% gather(Compound, intensity, -ID, na.rm = T)
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



