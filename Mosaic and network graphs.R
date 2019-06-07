#subset human, rat and mouse data only from interventions
interventions <- read.csv("data/interventions.csv")
inter.humanrat <- droplevels(subset(interventions, Species %in% c("Human", "Rat")))
inter1 <- interventions %>% filter(Species == "Human" | Species == "Rat")
                       
library(vcd)
library(labeling)

#plot of study designs
mosaic(~ Species + biofluids + Dose.type, data=inter.humanrat, highlighting="Species")

#plot of subclasses studied (how to change text size?)
mosaic(~ Subclass.studied + Species, data=inter.humanrat, highlighting="Species",
       direction=c("h", "v"), keep_aspect_ratio = F,
              labeling = labeling_border(rot_labels = c(0,0,0,0) ),
       labeling_args=list(gp_labels=gpar(fontsize=10), 
                          gp_varnames=gpar(fontsize=12)))

#reorder factor levels
inter.humanrat$Subclass.studied<-factor(inter.humanrat$Subclass.studied,
        levels=c("Flavonols", "Flavanones", "Anthocyanins", "Isoflavones", "Flavanols",
                 "HCA", "HBA", "Lignans", "Tyrosols", "Stilbenes"))

#with base graphics
plot(inter.humanrat$Subclass.studied, inter.humanrat$Species, 
     xlab=("Subclass studied"), ylab=("Species"), c = 0.8)

# Alluvial plot
# From https://matthewdharris.com/2017/11/11/a-brief-diversion-into-static-alluvial-sankey-diagrams-in-r/
dat_raw <- inter1 %>% select(Species, biofluids, Dose.type)
dat <- dat_raw %>%
        group_by(Species, biofluids, Dose.type) %>%
        summarise(freq = n()) %>%
        ungroup()

# Define colours
A_col <- "darkorchid1"
B_col <- "darkorange1"
C_col <- "skyblue1"
alpha <- 0.7 # transparency value
fct_levels <- c("A","C","B")

# Plot with alluvial
library(alluvial)
alluvial(dat[, 1:3], freq = dat$freq, gap.width=0.05, blocks = T,
         col= c(rep(A_col, 9), rep(B_col, 9)),
         border = c(rep(A_col, 9), rep(B_col, 9)),
         cw = 0.1, xw = 0.2, #alpha = 0.7,
         axis_labels = c("Species", "Biofluid", "Dose type")
         )

#-----------------------------------------------------------------------------------------

library(igraph)
# Network plot of Phenol-Explorer metabolism reactions
# Normal network graph
pathsnames <- read.csv("data/pathsnames.csv", header=T)
metadata <- read.csv("data/metadata.csv", header=T)

#create graph object from data set
g <-graph.data.frame(pathsnames,directed=TRUE)

#remove unnecessary margins
par(mar=c(1,1,1,1))

#plot graph
plot(g, 
     layout = layout.fruchterman.reingold, 
     #layout = layout.graphopt,
     vertex.size=3, 
     edge.arrow.size=0.2, 
     vertex.label.cex = 0.5,  # Slightly smaller font
     vertex.label.dist = 0.1,  # Offset the labels
     vertex.label.color = "black",
     vertex.label.family = "sans",
     vertex.color = "pink",
     vertex.frame.color="blue",
     vertex.label.font = 1)

# Circular network graph
#create graph object from data set
g <- graph.data.frame(pathsnames, vertices = metadata, directed= T )

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

