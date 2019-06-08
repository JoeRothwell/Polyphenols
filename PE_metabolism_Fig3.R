pkplasma <- read.csv("pk_plasma.csv", header=T)
interventions <- read.csv("interventions.csv", header=T)

#find median cmax: all doses
cmaxhuman <- pkplasma %>%
  filter(species == "Human") %>%
  left_join(interventions, by = "inter") %>%
  select(Dose.type, compound, cmax.um) %>%
  filter(!(is.na(cmax.um)))

#Food or pure compound only
cmaxhuman1 <- cmaxhuman %>% filter(Dose.type =="Food" | Dose.type =="Experimental food")
cmaxhuman2 <- cmaxhuman %>% filter(Dose.type =="Pure compound") # %>% arrange(desc(cmax.um))

median(cmaxhuman$cmax.um)
median(cmaxhuman1$cmax.um)

#----------------------------------------------------------------------------
#find median tmax for food and experimental food only
tmaxhuman <- pkplasma %>%
  filter(species == "Human") %>%
  left_join(interventions, by = "inter") %>%
  select(compound, tmax) %>%
  filter(Dose.type =="Food" | Dose.type =="Experimental food") %>%
  filter(!(is.na(cmax.um)))
#arrange(desc(cmax.um))

#----------------------------------------------------------------------------
#ecdf Cmax for figure

cmaxsource <- pkplasma %>%
  filter(species == "Human" | species == "Rat") %>%
  left_join(interventions, by = "inter") %>%
  mutate(purepp = ifelse(Dose.type == "Pure compound","PurePP", "Food"))

library(ggplot2)
library(scales)
source("Base breaks function.R")

fig3a <- ggplot(cmaxsource, aes(x=cmax.um, linetype=purepp)) + 
  stat_ecdf(aes(size=species)) + theme_bw(base_size=10) + 
  #scale_colour_manual(values=c("black", "darkgrey")) +
  #scale_linetype_manual(values=c("solid", "dashed")) +
  ##default line size (thickness) seems to be 0.5
  scale_size_manual(values=c(1,0.5)) +
  theme(axis.title.y=element_blank(), legend.position="none",
        panel.border=element_rect(colour="black", fill=NA),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
  scale_x_continuous(name=(bquote("C" [max]~ "("*mu*"M)")), trans=log_trans(), breaks=base_breaks(), 
                     labels=prettyNum) +
  scale_y_continuous(expand=c(0,0), limits=c(-0.02,1)) #+ 
#annotate("text", x=0.08, y=0.75, size=3, label="Human") + ggtitle("A")
#annotate("text", x=1, y=0.4, size=3, label="Rat") + 
#ggtitle("A")

ggsave("ECDF Cmax.png", width=100, height=65, unit="mm", plot=fig3a)

ggsave("ECDF Cmax.svg", width=100, height=65, unit="mm", plot=fig3a)

library(latticeExtra)
ecdfplot(~cmax.um, groups=species, data=cmaxsource, xlab= "Cmax (um)", ref=F, auto.key=F,
         scales=list(x=list(log=10, equispaced.log=F, tick.number=6)))

#-------------------------------------------------------------------------------------

#draws a frequency polygon for Tmax or Cmax by species

#subset human and rat data from pk_plasma
pkplasma.hr <- subset((read.csv("pk_plasma.csv")), species == "Human" | species == "Rat")

fig3b <- ggplot(pkplasma.hr, aes(x=tmax, size=species)) + geom_freqpoly() +
  scale_size_manual(values=c(1,0.5)) +
  theme_bw(base_size=10) + 
  theme(legend.position="none", panel.border=element_rect(colour="black", fill=NA),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
  annotate("text", x=2, y=48, size=3, label="Rat") + 
  annotate("text", x=4, y=39, size=3, label="Human") + 
  scale_x_continuous(name=(bquote("T"[max]~ "(h)")), expand=c(0.03,0.03), breaks=c(0,5,10,15,20,25)) +
  scale_y_continuous(name="No. values", expand=c(0.05,0.05)) #+
#ggtitle("B")

ggsave("frequency polygon tmax.png", width=100, height=65, unit="mm", plot=fig3b)

ggsave("frequency polygon tmax.svg", width=100, height=65, unit="mm", plot=fig3b)
