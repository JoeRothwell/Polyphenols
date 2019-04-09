#Food group by analysis type

library(ggplot2)
library(dplyr)
rfs <- read.csv("raw_data.csv", header=T)
aggregation2 <- rfs %>% group_by(food.group, analysis.type) %>% summarise(Count=n())

supp1a <- ggplot(aggregation2, aes(x=food.group, y=Count, fill=analysis.type)) + 
  geom_bar(stat="identity", position="dodge", colour="black", width=0.75) +
  theme_minimal() +
  theme(legend.position=c(0.8,0.8), axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_fill_brewer(palette="Greys") +
  scale_y_continuous(name="Number of retention factors", expand=c(0,0),
                     breaks=c(0,200,400,600,800,1000,1200)) +
  scale_x_discrete(limits=c("Vegetables", "Fruits", "Seeds", "Non-alc bev", "Cereals", "Alc bev", "Oils"),
                   labels=c("Vegetables", "Fruits", "Seeds", "Non-\nalcoholic\nbeverages",
                            "Cereals", "Alcoholic\nbeverages", "Oils"))

ggsave("Barchart agg vs raw RFs.png", width=133, height=100, unit="mm", plot=supp1a)

#Food group by compound class
library(plyr)
cmpdclasses <- ddply(rfs, c("food.group", "compound.class"), summarise, countrf=length(mean.rf), .drop=FALSE)

#dplyr: no drop=F at the moment!
cmpdclasses2 <- rfs %>% group_by(food.group, compound.class) %>% summarise(countrf=n())

#reorder factor levels
cmpdclasses$compound.class <- factor(cmpdclasses$compound.class, 
  levels=c("Flavonoids", "Phenolic acids", "Total phenolics", "Other polyphenols"))

#food group order
grps <- c("Vegetables", "Fruits", "Seeds", "Non-alc bev", "Cereals", "Alc bev", "Oils")

library(grid)

supp2 <- ggplot(cmpdclasses2, aes(x=food.group, y=countrf, fill=compound.class)) + 
  geom_bar(stat="identity", position="dodge", colour="black", width=0.7) + theme_minimal() +
  theme(legend.position=c(0.75,0.8)) + scale_fill_brewer(palette="Greys") +
  guides(fill=guide_legend(title=NULL)) +
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), legend.key.height=unit(4, "mm")) +
  scale_y_continuous(name="Aggregated retention factors", expand = c(0,0)) +
  scale_x_discrete(limits=grps, labels=c("Vegetables", "Fruits", "Seeds", "Non-\nalcoholic\nbeverages",
                            "Cereals", "Alcoholic\nbeverages", "Oils"))

ggsave("bar RFs by food group.png", width=133, height=100, unit="mm", plot=supp2)

#with lattice
library(lattice)
my.plot <- barchart(reorder(food.group,countrf) ~ countrf, groups=compound.class,
                    auto.key=T, data=cmpdclasses2, origin=0, xlab="No. retention factors")
trellis.device(device="png", filename="RFs bar lattice.png")
print(my.plot)
dev.off()

#Food group by compound class 2
library(reshape2)
foodgp <- dcast(rfs, food.group + compound.class ~ ., drop=F)

ggplot(foodgp, aes(x=food.group, y=., fill=compound.class)) + 
  geom_bar(stat="identity", position="dodge", colour="black", width=0.7) + theme_minimal() +
  theme(legend.position=c(0.75,0.8)) + scale_fill_brewer(palette="Greys") +
  guides(fill=guide_legend(title=NULL)) +
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), legend.key.height=unit(4, "mm")) +
  scale_y_continuous(name="Aggregated retention factors", expand = c(0,0)) +
  scale_x_discrete(limits=grps, labels=c("Vegetables", "Fruits", "Seeds", "Non-\nalcoholic\nbeverages",
                            "Cereals", "Alcoholic\nbeverages", "Oils"))
ggsave("food group bar.png")

#Bar food groups

#summarise rfs, not dropping the zero categories
rfsbreakdown <- ddply(rfs, c("food.group", "process.type"), summarise, rfcount=length(mean.rf), .drop=FALSE)
#dplyr doesn't drop empty levels
rfsbreakdown <- rfs %>% group_by(food.group, process.type) %>% summarise(rfcount=n())

#beside
fig2 <- ggplot(rfsbreakdown, aes(x=reorder(food.group,rfcount), y=rfcount, fill=process.type)) + 
  geom_bar(stat="identity", position="dodge", colour="black", width=0.8, show_guide=FALSE) +
  theme_minimal(base_size=9) + guides(fill=guide_legend(title=NULL)) + 
  theme(legend.position=c(0.75,0.9), legend.justification=c(0.5,0.5), axis.ticks.x=element_blank(),
        axis.title.x=element_blank(), axis.text.x=element_text(), panel.grid.major=element_blank()) +
  scale_fill_brewer(palette="Greys") +
  labs(x="Food group", y="No. retention factors") + theme(legend.key.height=unit(3, "mm")) +
  scale_y_continuous(expand = c(0,0), limits=c(0,350)) +
  scale_x_discrete(limits=grps, "food.group", labels=c("Vegetables", "Fruits", "Seeds", 
                  "Non-\nalcoholic\nbeverages", "Cereals", "Alcoholic\nbeverages", "Oils"))

ggsave("fig2revised.pdf", width=80, height=70, unit="mm", dpi=600)

#stacked
fig2 <- ggplot(rfsbreakdown, aes(x=food.group, y=rfcount, fill=process.type)) + 
  geom_bar(stat="identity", colour="black", width=0.6) + theme_classic() + 
  guides(fill=guide_legend(title=NULL)) + 
  theme(legend.position=c(1,1), legend.justification=c(1,1)) +
  scale_fill_brewer(palette="Greys") +
  labs(x="Food group", y="No. retention factors") +
  theme(axis.text.x=element_text(angle=15, hjust=0.5), panel.border=element_rect(fill=NA), legend.key.height=unit(4, "mm")) +
  scale_y_continuous(expand = c(0,0), limits=c(0,500)) +
  scale_x_discrete(limits=grps)

ggsave("fig2.png", width=120, height=90, unit="mm", dpi=600, plot=fig2)

#Bar RFs by food group (base graphics)
barplot(table(rfs$process.type, rfs$food.group)[, c(7,3,6,4,2,1,5)], beside=T)

#Bar RF by process (ordered)
bar <- rfs %>% group_by(process.full) %>% summarise(ct=n()) %>%
       ggplot(aes(x=reorder(process.full, ct), y=ct)) + geom_bar(stat="identity", fill="grey") + 
       theme_bw() + coord_flip() + labs(x="Process", y="Aggregated retention factors")

ggsave("process barchart ordered.png", width=6, height=6, unit="in", plot=bar)

#--------------------------------------------------------------------------------------------------------

#dot plot processes
#to summarise RF data by process with count of raw and aggregated RFs and median and IQR of aggregated RFs
process.summary1 <- ddply(rfs, c("process.grouped", "process.type"), summarize, countraw=sum(nraw.rf),
    countagg=length(rf.id), medianagg=round(median(mean.rf), 2), iqragg=round(IQR(mean.rf), 2))

#to reorder by volume of data (count of raw RFs)
process.summary <- process.summary[order(process.summary$countraw, decreasing=TRUE), ] 

fig1 <- ggplot(process.summary1, aes(x=countagg, y=reorder(process.grouped, countagg))) + 
  geom_segment(aes(yend=process.grouped), xend=0, colour="black") + 
  geom_point(size=2.5, aes(shape=process.type)) + theme_bw(base_size=10) +
  scale_x_continuous(expand=c(0,0), limits=c(0,255)) +
  xlab("No. retention factors") +
  theme(axis.title.y=element_blank(), legend.position="none",
        panel.border=element_rect(colour="black"), axis.ticks.y=element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour="black", linetype="dotted"),
        strip.background = element_rect(colour="black", size=0.5)) +
  facet_grid(process.type ~ ., scales="free_y", space="free_y")

ggsave("fig1 revised.pdf", width=80, height=100, unit="mm", dpi=600, plot=fig1)

#----------------------------------------------------------------------------------------------------

#ECDF RFs by method
#remove pH differential method by subsetting
rfs <- read.csv("raw_data.csv", header=T)
rfs.nophdiff <- subset(rfs, !(method %in% "pH diff."))

ggplot(rfs.nophdiff, aes(x=mean.rf, colour=method)) + theme_bw() + stat_ecdf() +
  theme(legend.position=c(0.8,0.5), axis.title.y=element_blank()) +
  scale_x_log10(name="Aggregated retention factor", breaks=c(0.01,0.05,0.1,0.5,1,5,10,50)) 

ggsave("Cumulative density by method.png", width=133, height=100, unit="mm")

library(latticeExtra)
#to exclude ph diff using logical vector rfs$method!="pH diff."

ecdfplot(~mean.rf, groups=method, data=rfs, scales=list(x=list(log=10, equispaced.log=F)),
         xlab="Mean retention factor", auto.key=T, subset=rfs$method!="pH diff.")