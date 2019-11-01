#Correlation heatmap for urine

#read in, rows as obs and columns as variables
urinecor <- read.csv("data/urinary PPs names.csv", header=T)
urinecor1 <- as.matrix(urinecor)
#check variances?
apply(urinecor1, 2, var)

#generate correlation matrix
tcor <- cor(urinecor1, use="pairwise.complete.obs", method="pearson")

library(corrplot)
corrplot(tcor, mar=c(1,1,1,1), is.corr=FALSE, method="color",
         shade.col=NA, tl.col="black", tl.srt=90, order="hclust", type="full",
         tl.cex=0.8, cl.cex=0.7)

#cluster by distance matrix
hc <- hclust(dist(tcor))

#plot changing margins: c(bottom, left, top, right)
par(mar=c(0, 4, 2, 2))
pp <- plot(hc, xlab="", sub="", main="")

#heatmap of urinary PPs v food intake
#read data in and convert to matrix
library(dplyr) #marked foods to filter out with _rm
library(readr)
food <- read_csv("data/Correlation PPs with foods abbrev.csv") %>%
  filter(Polyphenol != "Gallocatechin") %>% select(-ends_with("rm"), -Polyphenol)

#draw heatmap in gplots exported at 7x5 inch:
library(gplots)
png(filename="Heatmap PPs food intake.png", width=130, height=200, res=200, units="mm")
heatmap.2(data.matrix(food), dendrogram="both", key=F, col=redblue(256), trace="none", 
          offsetRow = 0.2, offsetCol = 0.2, margins=c(6,6))
dev.off()

#read data in and convert to matrix
food <- read.csv("data/Correlation PPs with foods.csv", row.names=1)
colfoods <- read.csv("data/food names heatmap.csv", header=F)

png(filename="Heatmap PPs food intake key.png", width=200, height=200, res=200, units="mm")
heatmap.2(data.matrix(food), dendrogram="both", key=T, col=redblue(256), trace="none", 
          offsetRow = 0.1, offsetCol = 0.1, margins=c(8,12), labCol = unlist(colfoods))
dev.off()

#heatmap in ggplot2
#set up a coloring scheme using colorRampPalette (from RBloggers post)
red=rgb(1,0,0); green=rgb(0,1,0); blue=rgb(0,0,1); white=rgb(1,1,1)
RtoWrange <- colorRampPalette(c(red, white ))
WtoGrange <- colorRampPalette(c(white, green)) 

library(ggplot2)
library(reshape2) 
melted_food <- melt(food)
ggplot(data = food, aes(x=reorder(Var1, value, FUN=mean),
      y=reorder(Var2, value, FUN=mean), fill=value)) + geom_tile() +
  theme(axis.text.x=element_text(angle=40, hjust=1)) + coord_flip()
  #see below for definitions of colour ranges
  scale_fill_gradient2(low=RtoWrange(100), mid=WtoGrange(100), high="gray")
ggsave("food intake heatmap ggplot.png", width=200, height=150, units="mm")

#----------------------------------------------------------------------------------------------

#Node and edge tables for cytoscape correlation map

ue <- read.csv("data/Urinary polyphenols_EPIC.csv")

nodetable <- ue %>% select(uPPc_1:uPPc_39) %>% gather(polyphenol, conc) %>%
  filter(!is.na(conc)) %>% 
  group_by(polyphenol) %>% 
  summarise(median(conc)) %>%
  inner_join(labs, by="polyphenol")

mat <- log(select(ue, uPPc_1:uPPc_39))
cmat <- cor(mat, use="pairwise.complete.obs", method="pearson")

edgetable <- melt(cmat, value.name="pcor")

write.csv(nodetable, file="node table.csv")
write.csv(edgetable, file="edge table.csv")

#-------------------------------------------------------------------------------------

#old: remove gallocatechin and some redundant foods and convert to matrix with base R
v <- colnames(food)
w <- which(v=="Wine" | v=="Cereal_products" | v=="Stalk_veg.sprouts" | v=="Leafy_veg" |
             v=="Legumes" | v=="Nuts.seeds" | v=="Fruiting_veg" | v=="Root_veg" |
             v=="Grain.pod_veg" | v=="Cabbages" | v=="Breakfast_cereals")
#new matrix: indices for excluded foods are 14,1,8,2,9,10,3,4,6,5,12
food1 <- as.matrix(food[ -(10), -(w)])
