#read in food and metabolite data from csv
bc <- read.csv("all_bc.csv", header=T)
fcmpds <- read.csv("food_cmpds.csv", header=T)

#subset data tables
bch <- subset(bc, species=="human")
bcr <- subset(bc, species=="rat")
bct <- subset(bc, hydroylsis=="no")
bchy <- subset(bc, hydroylsis=="yes")
bcp <- subset(bc, biofluid=="plasma")
bcu <- subset(bc, biofluid=="urine")

#make vectors of the CICON for Venn function
hmetabs <- bch$cicon2
rmetabs <- bcr$cicon2
trmetabs <- bct$cicon2
hymetabs <- bchy$cicon2
allfcmpds <- fcmpds$CICON2
allmetabs <- bc$cicon2
pmetabs <- bcp$cicon2
umetabs <- bcu$cicon2

#draws a Venn diagram and exports image file. Give two vectors with labels for each.
library(VennDiagram)
#Venn 1: Human vs rat metabolites
venn.diagram(list(H=hmetabs, R=rmetabs), "Venn1.svg", imagetype = "svg", #fill = "grey",
             height=2.5, width=2.5, units="mm", main="A) Metabolites by species", main.pos=c(0.5, 1))

#Venn 2: Urine vs plasma metabolites
venn.diagram(list(U=umetabs, P=pmetabs), "Venn2.svg", imagetype = "svg", #fill = "grey",
             height=2.5, width=2.5, units="mm", main="B) Metabolites by biofluid", main.pos=c(0.5, 1))

#Venn 3: Hydrolysed vs non-hydrolysed metabolites
venn.diagram(list(N=trmetabs, Y=hymetabs), "Venn3.svg", imagetype = "svg", #fill = "grey",
             height=2.5, width=2.5, units="mm", main="C) Metabolites by hydrolysis", main.pos=c(0.5, 1))

#Venn 4: Food compounds vs true metabolites
venn.diagram(list(FC=allfcmpds, M=trmetabs), "Venn4.svg", imagetype = "svg", #fill = "grey",
             height=2.5, width=2.5, units="mm", main="D) All polyphenols by type", main.pos=c(0.5, 1))

#Publications

#extracts the cipub IDs as vectors from the interventions data frame
interventions <- read.csv("interventions.csv", header=T)

humanpubs <- (subset(interventions, species=="human"))$cipub
ratpubs <- (subset(interventions, species=="rat"))$cipub
mousepubs <- (subset(interventions, species=="mouse"))$cipub

#draws a Venn diagram and exports it as a tiff. Give two vectors with labels for each.
venn.diagram(list(H=humanpubs, R=ratpubs, M=mousepubs), "publication species venn3.tiff",
             height=90, width=90, units="mm") 

