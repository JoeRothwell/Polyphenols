#subset human, rat and mouse data only from interventions
inter.humanrat <- droplevels(subset(interventions, Species %in% c("Human", "Rat")))
                       
library(vcd)
library(labeling)

#plot of study designs
mosaic(~ Species + Biofluids + Dose.type, data=inter.humanrat, highlighting="Species")

#plot of subclasses studied (how to change text size?)
mosaic(~ Subclass.studied + Species, data=inter.humanrat, highlighting="Species",
       direction=c("h", "v"), keep_aspect_ratio=FALSE)

#reorder factor levels
inter.humanrat$Subclass.studied<-factor(inter.humanrat$Subclass.studied,
        levels=c("Flavonols", "Flavanones", "Anthocyanins", "Isoflavones", "Flavanols",
                 "HCA", "HBA", "Lignans", "Tyrosols", "Stilbenes"))

#with base graphics
plot(inter.humanrat$Subclass.studied, inter.humanrat$Species, 
     xlab=("Subclass studied"), ylab=("Species"))

#-----------------------------------------------------------------------------------------

library(igraph)
#create graph object from data set
g <-graph.data.frame(pathsnames,directed=TRUE)

#remove unnecessary margins
par(mar=c(1,1,1,1))

#plot graph
plot(g,layout=layout.fruchterman.reingold, vertex.size=3, edge.arrow.size=0.5, 
     vertex.label.cex  = 0.4,  # Slightly smaller font
     vertex.label.dist  = 0.1,  # Offset the labels
     vertex.label.color = "black",
     vertex.label.family = "sans",
     vertex.color = "pink",
     vertex.frame.color="blue",
     vertex.label.font = 1)


#Circular network graph
pathsnames <- read.csv("pathsnames.csv", header=T)
metadata <- read.csv("metadata.csv", header=T)

#create graph object from data set
g <- graph.data.frame(pathsnames, vertices=metadata, directed=TRUE)

#remove unnecessary margins
par(mar=c(1,1,1,1))

V(g)$color <- ifelse(V(g)$Type == 1, "orange", "green")

#plot graph
plot(g, layout=layout.circle, vertex.size=3, edge.arrow.size=0.2,
     vertex.label.family = "sans",
     vertex.label.cex  = 0.4,
     vertex.label.dist  = 0.1,  # Offset the labels
     vertex.label.color = "grey",
     vertex.frame.color="black",
     vertex.label.font = 1)

