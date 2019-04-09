#Scatter PK graphs
library(dplyr)
source("Base breaks function.R")

pk.plasma <- read.csv("pk_plasma.csv", header=T)
pk <- pk.plasma %>% filter(species == "Human" | species == "Rat")
pk2 <- pk %>% group_by(subclass) %>% filter(n() > 10 | subclass == "Lignans")

abbrev <- c("Hydroxycinnamic acids" = "Hydroxy-\ncinnamic\nacids", "Isoflavonoids"="Isoflavones", 
            "Hydroxybenzoic acids" = "Hydroxy-\nbenzoic\nacids")

library(ggplot2)
library(grid)
library(scales)

#Tmax and thalf by subclass and species
#panel.margin increases space between facets
fig5c <- ggplot(pk2, aes(x=subclass, y=tmax, shape=species)) +
  geom_boxplot(fill="dodgerblue") + guides(shape=FALSE) + coord_flip() +
  geom_jitter(position=position_jitter(w=0.2, h=0.1), alpha=0.7) +
  xlab("Tmax (hours)") + scale_shape_manual(values=c(20, 17)) +
  facet_grid(. ~ species, as.table=FALSE) + theme_bw(base_size=10) +
  theme(axis.title.y=element_blank(), panel.margin = unit(1, "lines")) +
  scale_x_discrete(labels=abbrev, limits=levels(droplevels(pk2$subclass))) +
  scale_y_continuous(trans=log_trans(), breaks=base_breaks(), labels=prettyNum) + ggtitle("C")

ggsave("tmax subclass facet horiz.svg", width=145, height=90, units="mm", plot=fig5c)

#reorder factor levels to reverse human and rat
pk.plasma3$species <- factor(pk.plasma3$species, levels=c("dog", "rat", "mouse", "pig", "rabbit", "human"))

fig4 <- ggplot(pk2, aes(x=subclass, y=tmax, shape=species)) + 
  geom_boxplot(varwidth=TRUE, outlier.size=0, colour="grey60") + guides(shape=FALSE) +
  geom_jitter(position=position_jitter(w=0.2, h=0.1), alpha=0.8) +
  xlab("Subclass") + ylab("Tmax (hours)") + scale_shape_manual(values=c(17, 20)) +
  facet_grid(species ~ ., as.table=FALSE) + theme_bw(base_size=9) +
  scale_x_discrete("subclass", labels=abbrev, limits=levels(droplevels(pk2$subclass))) +
  scale_y_log10(breaks=c(0.1,0.5,1.0,5,10)) + theme(panel.border =element_rect(fill=NA))

ggsave("tmax subclass facet2.png", width=150, height=110, units="mm", plot=fig4)

fig3c <- ggplot(pk2, aes(x=subclass, y=tmax, shape=species)) + coord_flip() +
  geom_jitter(position=position_jitter(w=0.1, h=0.3)) + 
  guides(shape=FALSE, colour=FALSE) +
  scale_shape_manual(values=c(1,6)) + theme_bw(base_size=10) + #ggtitle("C") +
  theme(axis.title.y=element_blank(), axis.ticks.y=element_blank(),
        panel.border=element_rect(color="black", fill=NA), panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  scale_x_discrete(labels=abbrev) +
  scale_y_continuous(name=(bquote("T"[max]~ "(h)")), trans=log_trans(),
                     breaks=base_breaks(), labels=prettyNum)

ggsave("tmax subclass.png", width=100, height=80, units="mm", plot=fig3c)

ggsave("tmax subclass.svg", width=100, height=80, units="mm", plot=fig3c)

fig4 <- ggplot(pk, aes(x=species, y=tmax, colour=subclass, shape=species)) +
  geom_jitter(position=position_jitter(w=0.1, h=0.15)) + scale_shape_manual(values=c(6,19)) +
  #geom_hline(yintercept=1, linetype="dashed") + ylim(0,5) + 
  theme_bw(base_size=10) + labs(x="Species", y="Tmax (hours)") + 
  scale_x_discrete(limits=c("human", "rat")) + scale_y_log10()

ggsave("tmax subclass scatter.png", width=150, height=100, unit="mm", plot=fig4)

thalf <- ggplot(pk, aes(x=subclass, y=thalf, shape=species)) +
  geom_jitter(position=position_jitter(w=0.05, h=0)) + coord_flip() + theme_bw(base_size=10) +
  guides(shape=FALSE) + ylab("Half life (hours)") + scale_shape_manual(values=c(6, 16)) +
  theme(axis.title.y=element_blank()) + scale_x_discrete("subclass", labels=c(
    "Hydroxyphenylpropanoic acids"="Hydroxyphenyl-\npropanoic acids")) +
  scale_y_log10(breaks=c(0.1,1,0.5,1,5,10,50))

ggsave("thalf scatter2.png", width=133, height=100, units="mm", plot=thalf)

fig4 <- ggplot(pk, aes(x=species, y=thalf, colour=subclass)) +
  geom_jitter(position=position_jitter(w=0.05, h=0)) +
  theme_bw(base_size=10) + labs(x="Species", y="Thalf (hours)") + 
  scale_x_discrete(limits=c("Human", "Rat"))

ggsave("thalf scatter.png", height=120, width=140, units="mm", plot=fig4)

#---------------------------------------------------------------------------------------------------

#Scientific notation

#subscript in axis labels: square brackets
#superscripts in square brackets followed by ~
#normal text must be enclosed by ''
# *mu~ causes greek letter mu then a space
# *mu* causes greek letter mu without space
xlab(   bquote("C" [max]~ "("*mu~"M)")    )

#---------------------------------------------------------------------------------------------------

#scatter plot of Tmax vs Cmax for all species
pk <- ggplot(pk.plasma, aes(x=tmax, y=cmax.um, shape=species)) + geom_point() +
  theme_bw() + scale_shape_manual(values=c(21,16,5,8,7,6)) + xlab("Tmax (h)") + ylab("Cmax (μM)")

pk + scale_y_log10()
ggsave("log species pk comparison2.png", height=100, width=150, unit="mm")

pk + scale_y_continuous(expand=c(0.02,0.02), limits=c(0,25)) + scale_x_continuous(expand=c(0.02,0.02), 
                                                                                  limits=c(0,12.5))
ggsave("species pk comparison ylim25.png", height=70, width=110, unit="mm")

#base graphics
plot(pk.plasma$tmax, log(pk.plasma$cmax.um))

#human only
pk <- ggplot(pk.plasma.human, aes(x=tmax, y=cmax.um, colour=subclass)) + geom_point() +
  xlim(0,12) + ylim(0,4) + theme_bw() + xlab("Tmax (hours)") + ylab("Cmax (μM)")

ggsave("human pk subclass.png", height=4, width=8, unit="in")

#---------------------------------------------------------------------------------------------------

#Urinary excretion (for supp data)
pkurine <- read.csv("pk_urine.csv", header=T)
#get hydrolysed data only for urinary recovery
ue <- subset(pkurine, !(hydrolysis %in% "None"))

suppdata <- ggplot(ue, aes(x=urinary.ex, y=reorder(subclass, urinary.ex), shape=species)) +
  geom_jitter(position=position_jitter(w=0.05, h=0.02)) + theme_bw(base_size=12) +
  guides(shape=FALSE) + 
  scale_shape_manual(values=c(1, 6)) +
  theme(axis.title.y=element_blank(), panel.grid.major=element_blank()) +
  scale_x_continuous(name="Percentage dose recovered") +
  scale_y_discrete("subclass", 
                   labels=c("Hydroxyphenylpropanoic acids"="Hydroxyphenyl-\npropanoic acids"))

ggsave("urinary recovery.png", width=150, height=90, units="mm", plot=suppdata)

#median and IQR of urinary recovery
median(ue$urinary.ex)
IQR(ue$urinary.ex)

#by species
pkurinehydrol <- subset(read.csv("pk_urine.csv"), !(hydrolysis %in% "None"))

fig4 <- ggplot(pkurinehydrol, aes(x=species, y=urinary.ex, colour=subclass)) +
  geom_jitter(position=position_jitter(w=0.05, h=0)) +
  theme_bw(base_size=10) + labs(x="Species", y="Urinary recovery %") + 
  scale_x_discrete(limits=c("human", "rat"))

ggsave("urinary excretion scatter.png", height=100, width=110, units="mm", plot=fig4)

#-------------------------------------------------------------------------------------------------------

#range plots (abandonned)

tmaxrange <- ggplot(pk.plasma.summary, aes(x=reorder(subclass, meantmax), y=meantmax, colour=species)) + 
  geom_point(position=position_dodge(width=0.5, height=0), size=2) + 
  geom_errorbar(aes(ymax=meantmax+sdtmax, ymin=meantmax-sdtmax),
                position=position_dodge(width=0.5, height=0), width=0.5) + 
  theme_bw() + coord_flip()

ggsave("tmax range plot.png", height=3, width=6, unit="in", plot=tmaxrange)

#second method
tmax.summary <- ddply(pk.plasma, c("subclass", "species"), summarise, mintmax=min(tmax), maxtmax=max(tmax))

tmaxrange <- ggplot(tmax.summary, aes(x=subclass, y=mintmax, colour=species)) + 
  geom_point(position=position_dodge(width=0.5, height=0), size=2) + theme_bw() + coord_flip() +
  geom_errorbar(aes(ymax=maxtmax, ymin=mintmax), position=position_dodge(width=0.5, height=0), width=0.5)

ggsave("tmax range plot.png", height=3, width=6, unit="in", plot=tmaxrange)
