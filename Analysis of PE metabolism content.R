#count detections of each metabolite by biofluid
library(dplyr)
all.bc <- read.csv("all_bc.csv", header=T)
detections <- all.bc %>% count(compound, biofluid) %>% ungroup() %>% arrange(desc(n))

#All combinations of dose type and species are listed with the number of interventions for each.
interventions <- read.csv("interventions.csv")
interventions %>% count(Dose.type, Species) %>% ungroup() %>% arrange(Dose.type, desc(n))

#Histogram of intervention duration
library(ggplot2)

fig1 <- ggplot(interventions, aes(x=duration.cat, fill=multiple.dose)) + 
  geom_histogram(position="dodge", colour="black", width=0.6) + theme_minimal() +
  theme(legend.position=c(0.75,0.7), axis.title.x=element_blank(), axis.text.x=element_text(),
        axis.ticks.x=element_blank()) + ylab("No. interventions") +
  scale_fill_brewer(palette="Greys", labels=c("Single dose", "Repeated dose"), guide = guide_legend(title=NULL)) +
  scale_x_discrete("duration.cat", labels=c("Less than 12 h" = "Less than\n12 hours",
        "12-24 h" = "12-24\nhours", "2-7 days" = "2-7\ndays", "More than 1 month" = "More than\n1 month"),
        limits=c("Less than 12 h", "12-24 h", "2-7 days", "More than 1 month")) +
  scale_y_continuous(expand=c(0,0))

ggsave("intervention length histogram.png", height=75, width=100, units="mm", plot=fig1)

#------------------------------------------------------------------------------------------------------
#dotplot of polyphenol sublass studied
library(plyr)
#to redo the below with dplyr
subclasses <- ddply(interventions, c("Species2", "Subclass.studied"),
                    summarise, Interventions=length(inter))

subclasses2 <- subset(subclasses, (Species %in% c("Human", "Rat")))

fig1 <- ggplot(subclasses.studied2, aes(x=Interventions, y=reorder(Subclass.studied, Interventions))) + 
  geom_segment(aes(yend=Subclass.studied), xend=0, colour="grey50") + 
  geom_point(size=2.5, aes(shape=Species)) + theme_bw() +
  scale_x_continuous(expand=c(0,0), limits=c(0,45)) + xlab("No. interventions") +
  scale_y_discrete(limits=c("Stilbenes", "Tyrosols", "Lignans", "HCA", "HBA", "Isoflavones",
                            "Flavonols", "Flavanones", "Flavanols", "Anthocyanins")) +
  theme(axis.title.y=element_blank(), legend.position="none",
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour="grey60", linetype="dashed")) +
  facet_grid(. ~ Species, scales="free_y", space="free_y")

ggsave("dot plot subclass studied horiz.png", width=100, height=67, unit="mm", plot=fig1)

#-----------------------------------------------------------------------------------------------

#Description of metabolite types

#subset true metabolites only
true.metabolites <- subset(true.metabolites, (hydrolysis %in% "no"))

#summarise data
summary.truemetabs <- ddply(true.metabolites, c("class", "class.no", "subclass", "subclass.count", "type"),
                            summarize, count=length(type))

#choose polyphenols for limits
PPs <- c("Anthocyanins", "Isoflavonoids", "Flavanols", "Flavanones", "Flavonols",
         "Dihydroflavonols", "Dihydrochalcones", "Chalcones", "Hydroxycinnamic acids", 
         "Hydroxybenzoic acids", "Hydroxyphenylpropanoic acids", "Hydroxyphenylacetic acids", 
         "Hydroxyphenylpentanoic acids", "Hydroxybenzoketones", "Hydroxybenzaldehydes", 
         "Hydroxyphenylalcohols", "Lignans", "Stilbenes", "Tyrosols", "Hydroxycoumarins", 
         "Other polyphenols", "Non-phenolic metabolites")

metabtypes <- ggplot(summary.truemetabs, aes(x=subclass, y=count, fill=type)) +
  geom_bar(stat="identity", width=0.7) + geom_bar(colour="black", show_guide=FALSE) +
  theme_bw() + coord_flip() + scale_fill_brewer(palette="Greys") +
  scale_y_continuous(expand = c(0,0), limits=c(0,60), breaks=c(0,10,20,30,40,50,60)) +
  labs(x="", y="Number of metabolites") + theme(legend.position=c(0.7,0.3))
guides(fill=guide_legend(title=NULL)) +
  scale_x_discrete(limits=PPs))

ggsave("true metabolite types.png", width=120, height=100, unit="mm", plot=metabtypes)

fig2 <- ggplot(true.metabolites, aes(x=subclass, fill=type)) +
  geom_bar(width=0.7, colour="black") + theme_bw() +
  scale_fill_brewer(palette="Greys") + coord_flip() +
  scale_y_continuous(expand = c(0,0), limits=c(0,60), breaks=c(0,10,20,30,40,50,60)) +
  labs(x="", y="Number of metabolites") + theme(legend.position=c(0.7,0.3))
guides(fill=guide_legend(title=NULL)) +
  scale_x_discrete(limits=PPs)

ggsave("metabolite types by subclass.png", width=150, height=100, units="mm", plot=fig2)

#-----------------------------------------------------------------------------------------

#base graphics: breakdown of metabolites by subclass and type. Need a matrix not df
library(reshape2)
types <- acast(subset(all.metabolites, hydrolysis == "no"), subclass ~ type)
#transpose and plot as horizontal bar chart
barplot(t(types), horiz=T, cex.names=0.5, las=1, legend=T)

#-------------------------------------------------------------------------------------------

#Bar graph metabolite classes
#summarize data
metabs <- ddply(all.metabolites, c("class", "type"), summarize, count=length(type), .drop=FALSE)

#reorder factor levels
metabs$type <- factor(metabs$type, levels=c("Aglycone", "Glucuronide", "Sulfate", "Glycoside"))
library(grid)

fig2 <- ggplot(metabs, aes(x=reorder(class, count), y=count, fill=type)) +
  geom_bar(stat="identity", position="dodge", width=0.8, colour="black") +
  theme_minimal() + 
  scale_x_discrete("class", labels=c("Flavonoids"="Flavonoids", "Phenolic acids"="Phenolic\nacids",
                                     "Other polyphenols" = "Other\npolyphenols", "Stilbenes"="Stilbenes",
                                     "Non-phenolic metabolites" = "Non-phenolic\nmetabolites"),
                            limits=c("Flavonoids", "Phenolic acids", "Lignans",
                                    "Stilbenes", "Other polyphenols")) +
  scale_fill_brewer(palette="Greys") +
  theme(legend.position=c(0.8,0.85), axis.title.x=element_blank(),
        axis.text.x=element_text(), axis.ticks.x=element_blank(), legend.key.height=unit(4, "mm")) +
  scale_y_continuous(expand = c(0,0), limits=c(0,100)) +
  ylab("Number of metabolites") +
  guides(fill=guide_legend(title=NULL))

ggsave("metabolite types beside.png", width=133, height=100, unit="mm", plot=fig2)

