#read in data and load packages
ue <- read.csv("Urinary polyphenols_EPIC.csv")
labs <- read.csv("uPPc codes.csv", header=T)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(grid)
source("Base breaks function.R")

#------------------------------------------------------------------------------------------

#subset data and make boxplot of concentrations
uemelt <- ue %>% select(Center, uPPc_1:uPPc_39) %>% gather(polyphenol, conc, -Center) %>% 
  filter(!is.na(conc)) %>% inner_join(labs, by="polyphenol")

#stat summary line gets median value for each epic centre and superimposes a point over the box
ggplot(uemelt, aes(x=reorder(ppname, conc, FUN=median), y=conc)) + 
  geom_boxplot(outlier.colour = "white", fill="grey", size=0.4) + 
  stat_summary(aes(shape=factor(Center)), fun.y=median, geom="point", size=1) + 
  scale_shape_manual(values=c(1,1,1,1,1,1,1,1,1)) +
  coord_flip() + theme_bw() +
  guides(shape=FALSE) +
  theme(axis.title.y=element_blank(), panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  scale_y_continuous(trans=log_trans(), breaks=base_breaks(), labels=prettyNum,
        name=(bquote("Median urinary excretion ("*mu*"g/24h)")))

ggsave("Boxplot PPs manuscript2.png", width=140, height=130, units="mm")
ggsave("Boxplot PPs manuscript.pdf", width=140, height=130, units="mm")

#------------------------------------------------------------------------------------------

#effects plot by country in black and white
ue %>% group_by(Country) %>% select(-(Idepic:RE_ENERGY)) %>%
  summarise_each(funs(median(., na.rm=T))) %>% gather(polyphenol, Excretion, -Country) %>%
  inner_join(labs, by="polyphenol") %>%
  ggplot(aes(x=reorder(Country, Excretion), y=Excretion, shape=class)) + 
  geom_line(aes(group=ppname)) + geom_point(aes(group=ppname)) +
  scale_y_continuous(name = "Median urinary excretion (ug/24h)", trans=log_trans(), 
                     breaks=base_breaks(), labels=prettyNum) + 
  xlab("Country") + theme_classic() +
  annotate("text", x=4.3, y=179, size=3, label="1") + # 4-Hydroxyphenylacetic acid
  annotate("text", x=4.3, y=65, size=3, label="2") + # Ferulic acid
  annotate("text", x=4.3, y=50, size=3, label="3") + # Vanillic acid
  annotate("text", x=4.3, y=40, size=3, label="4") + # 3-Hydroxyphenylacetic acid
  annotate("text", x=4.3, y=30, size=3, label="5") + # Homovanillic acid
  annotate("text", x=4.3, y=25, size=3, label="6") + # 4-Hydroxybenzoic acid
  annotate("text", x=4.3, y=20, size=3, label="7") + # 3,5-Dihydroxyphenylpropionic acid
  annotate("text", x=4.3, y=13, size=3, label="8") + #
  annotate("text", x=4.3, y=8,  size=3, label="9") + #
  annotate("text", x=4.5, y=7,  size=3, label="10") + #
  annotate("text", x=4.3, y=5.5, size=3, label="11") +
  annotate("text", x=4.5, y=5,   size=3, label="12") +
  annotate("text", x=4.3, y=4.5, size=3, label="13") +
  annotate("text", x=4.5, y=4,   size=3, label="14") +
  annotate("text", x=4.3, y=3,   size=3, label="15") +
  annotate("text", x=4.3, y=2.5, size=3, label="16") +
  annotate("text", x=4.5, y=2.2, size=3, label="17") +
  annotate("text", x=4.3, y=1.7, size=3, label="18") +
  annotate("text", x=4.5, y=1.5, size=3, label="19") +
  annotate("text", x=4.3, y=1.2, size=3, label="20") +
  annotate("text", x=4.3, y=0.8, size=3, label="21") +
  annotate("text", x=4.3, y=0.6, size=3, label="22") +
  annotate("text", x=4.5, y=0.55, size=3, label="23") +
  annotate("text", x=4.3, y=0.5,  size=3, label="24") +
  annotate("text", x=4.5, y=0.45, size=3, label="25") +
  annotate("text", x=4.3, y=0.4,  size=3, label="26") +
  annotate("text", x=4.5, y=0.36, size=3, label="27") +
  annotate("text", x=4.3, y=0.27, size=3, label="28") +
  annotate("text", x=4.3, y=0.16, size=3, label="29") +
  annotate("text", x=4.3, y=0.13, size=3, label="30") +
  annotate("text", x=4.3, y=0.11, size=3, label="31") +
  annotate("text", x=4.3, y=0.08, size=3, label="32") +
  annotate("text", x=4.3, y=0.065, size=3, label="33") +
  annotate("text", x=4.3, y=0.05, size=3, label="34") +
  theme(legend.position="none", axis.title.x=element_blank())
  #theme(legend.key.height=unit(0.3, "cm"))
  
ggsave("Line graph PP excretion.png", height=150, width=110, units="mm")

#------------------------------------------------------------------------------------------

#correlation heatmap of metabolites. First make names
foodcor <- read.csv("Correlations PP with foods new.csv")
ppnames <- as.vector(foodcor$Polyphenol)

mat <- log(select(ue, uPPc_1:uPPc_39))
cmat <- cor(mat, use="pairwise.complete.obs", method="pearson")
rownames(cmat) <- ppnames
colnames(cmat) <- vector()

library(corrplot)
png(filename="Correlation heatmap PPs2.png", width=200, height=200, res=300, units="mm")
corrplot(cmat, is.corr=FALSE, method="color", tl.col="black", order="hclust", type="full",
         tl.cex=1, tl.pos="tl", cl.pos = "b")
dev.off()

#correlation heatmap of food intake and metabolites
mat <- data.matrix(foodcor[, 2:12])
foodnames <- c("Olive oil", "Red wine", "Coffee", "Tea", "Herbal tea",
               "Citrus juices", "Bread, non-white", "Apple & pear",
               "Citrus fruits", "Chocolate & candy", "Onion & garlic") #in order of original df
library(gplots)
png(filename="Heatmap PPs food intake new.png", width=130, height=220, res=300, units="mm")
heatmap.2(mat, dendrogram="both", key=F, col=redblue(256), trace="none", 
          offsetRow = 0, offsetCol = 0.1, margins=c(9,11), 
          labRow = ppnames, labCol = foodnames)
dev.off()

#------------------------------------------------------------------------------------------

#Effects plot in colour
library(readr)
profilespp <- read_csv("Urinary PPs names.csv", col_names=T)
#gather argument says "gather all variables except country calling the new key column "compound"
#and the new value column "excretion"
meds <- profilespp %>% group_by(country) %>% select(-(Idepic:bmi)) %>% 
  summarise_each(funs(median(., na.rm=T))) %>% gather(Compound, Excretion, -country)

source("Base breaks function.R")
ggplot(meds, aes(x=reorder(country, Excretion), y=Excretion, colour=Compound)) + 
  geom_line(aes(group=Compound)) +
  geom_point(aes(group=Compound)) + #ylab("Median urinary excretion (ug)") +
  scale_y_continuous(name = "Median urinary excretion (ug)", trans=log_trans(), breaks=base_breaks(),
                     labels=prettyNum) + 
  xlab("Country") + theme_classic() +
  theme(legend.key.height=unit(0.3, "cm"),
#added as theme_classic no longer draws axes
  axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))

ggsave("Line graph PP excretion.png", height=120, width=200, units="mm")



