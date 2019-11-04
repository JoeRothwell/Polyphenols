# Analysis of RF data. Originally 6 scripts, now concatenated

# Subsetting

#to subset a data frame
rfs <- read.csv("data/raw_data.csv", header=T)
rfs.chrom1 <- subset(rfs, method=="Chrom.") #or
rfs.chrom2 <- rfs[rfs$method=="Chrom.", ]

#to drop the unused factor levels
droplevels(rfs.chrom)

#combined
rfs.chrom <- droplevels(subset(rfs, method=="Chrom"))

#subset by index
rfs.filter<-rfs[rfs$process %in% c("baked", "blanched", "boiled", "grilled", "microw.", 
                                   "pr. boiled", "pr. steamed", "roasted"),]
rfs.filter$process<-factor(rfs.filter$process)
for (i in ncol(rfs.filter)){}
for (i in ncol(rfs.filter)){rfs.filter[,i]<-factor(rfs.filter[,i])}

#-------------------------------------------------------------------------------------------------

rfs.filter <- rfs[rfs$process %in% c("baked", "blanched", "boiled", "grilled", "microw.",
                         "pr. boiled", "pr. steamed", "roasted"),]

#using a loop (not sure if good idea)
rfs.filter$process <- factor(rfs.filter$process)
for (i in ncol(rfs.filter)){}
for (i in ncol(rfs.filter)){rfs.filter[, i] <- factor(rfs.filter[,i])}

#using subset function with ! to exclude storage and industrial processing
rfs.domcook <- subset(rfs, !(process.type %in% c("Storage", "Ind proc")))
                      
#same result
rfs.domcook <- subset(rfs, (process.type %in% "Dom cook"))
                      
#RFs with dplyr    
library(dplyr)
                      
#makes a tbl data frame for dplyr
rfs2 <- tbl_df(rfs)
                      
#subset only certain columns
rfs.reduced <- select(rfs2, food, process, compound, mean.rf)
                      
#subset only certain rows
rfs.reduced.boil <- filter(rfs.reduced, process == "Boiled")
                      
rfs.byprocess <- group_by(rfs2, process)
summarise(rfs.byprocess, median(mean.rf), IQR(mean.rf))

#------------------------------

# RF summary tables

#By compound
library(plyr)
library(dplyr)
rfs <- read.csv("raw_data.csv", header=T)

#subset defined and grouped compounds, excluding storage processes
cmpds.def <- rfs %>% filter(process.type != "Storage" & method == "Chrom.") %>% group_by(compound) %>% 
  summarise(nraw = sum(nraw.rf), nagg=length(rf.id), medagg=median(mean.rf), iqragg = IQR(mean.rf))

cmpds.grp <- rfs %>% filter(process.type != "Storage" & analysis.type == "Grouped") %>% group_by(compound) %>% 
  summarise(nraw = sum(nraw.rf), nagg=length(rf.id), medagg=median(mean.rf), iqragg = IQR(mean.rf))

#merge to give separate columns for defined and grouped, then order by sum of raw data count columns
cmpds.all <- merge(cmpds.def, cmpds.grp, by="compound", suffixes = c("_def", "_group"), 
                   all.x=TRUE, all.y=TRUE)
#with dplyr
cmpds.all2 <- full_join(cmpds.def, cmpds.grp, by="compound") %>% arrange(desc(nraw.x + nraw.y))

cmpds.all <- arrange(compounds.allpp, desc(countraw_def+countraw_group))

write.csv(cmpds.all, file="Stats RFs by compound2.csv")

#By food. summarise non-hydrolysed analyses only using plyr
food.def <- ddply(chrom, "food", summarize, countraw=sum(nraw.rf), 
                  countagg=length(rf.id), medianagg=round(median(mean.rf), 2), iqragg=round(IQR(mean.rf), 2))

#same for grouped analyses 
food.grp <- ddply(grp, "food", summarize, countraw=sum(nraw.rf), 
                  countagg=length(rf.id), medianagg=round(median(mean.rf), 2), iqragg=round(IQR(mean.rf), 2))

#merge to give separate columns for defined and grouped, then order by sum of raw data count columns
food.all <- merge(food.def, food.grp, by="food", suffixes = c("_def", "_group"), all.x=TRUE, all.y=TRUE)
food.all <- arrange(food.all, desc(countraw_def + countraw_group))

write.csv(foods.allpp, file="Stats foods defined and grouped RFs.csv")

#By process
#summarise non-hydrolysed analyses only and order
chr.pool <- ddply(rfs.chrom, "process.full", summarize, countraw=sum(nraw.rf), 
                  countagg=length(rf.id), medianagg=round(median(mean.rf), 2), iqragg=round(IQR(mean.rf), 2))
chr.pool <- chrom.pooled[order(chrom.pooled$countraw, decreasing=TRUE), ]

#summarise grouped analyses only and order
grp.pool <- ddply(rfs.grouped, "process.full", summarize, countraw=sum(nraw.rf),
                  countagg=length(rf.id), medianagg=round(median(mean.rf), 2), iqragg=round(IQR(mean.rf), 2))
grp.pool <- grp.pool[order(grp.pool$countraw, decreasing=TRUE), ]

#merge to give separate columns for defined and grouped, then order by sum of raw data count columns
all.pool <- merge(chr.pool, grp.pool, by="process.full",
                  suffixes = c("_chrom", "_group"), all.x=TRUE, all.y=TRUE)
all.pool <- all.pool[order((all.pool$countraw_chrom + all.pool$countraw_group), decreasing=TRUE), ]

#write to csv
write.csv(all.pool, file="Stats for defined and grouped RFs.csv")


# to summarise RF data by process with count of raw and aggregated RFs and median and IQR of aggregated RFs
procsum1 <- ddply(rfs, c("process.grouped", "process.type"), summarize, countraw=sum(nraw.rf),
                  countagg=length(rf.id), medianagg=round(median(mean.rf), 2), iqragg=round(IQR(mean.rf), 2))

#to reorder by volume of data (count of raw RFs)
procsum <- process.summary[order(process.summary$countraw, decreasing=TRUE), ] 

#merge summarised data frames (for processes)
procsum1 <- ddply(rfs.chrom, "process.full", summarize, countraw=sum(nraw.rf),
                  countagg=length(rf.id), medianagg=round(median(mean.rf), 2), iqragg=round(IQR(mean.rf), 2))

procsum2 <- ddply(rfs.chromhydrol, "process.full", summarize, countraw=sum(nraw.rf),
                  countagg=length(rf.id), medianagg=round(median(mean.rf), 2), iqragg=round(IQR(mean.rf), 2))

#merges data frames by a common column, adding a suffix to other columns from
#both tables and specifying to keep all rows where data are not present in both tables
split <- merge(procsum1, procsum2, by="process.full", suffixes = c("1", "2"), all.x=TRUE, all.y=TRUE)

#order new table by sum of countraw
split[order((split$countraw1 + split$countraw2), decreasing=TRUE), ]

#-------------------------------------------------------------------

# Scatter plots of RFs

#RFs for black bean by compound on a log scale
rfs <- read.csv("data/raw_data.csv", header=T)
bb <- read.csv("data/black_bean2.csv", header=T)
library(ggplot2)
library(grid)
library(scales)
source("Base breaks function.R")
theme_set(theme_bw(10))

fig6a <- ggplot(bb, aes(x=mean.rf, y=reorder(compound, mean.rf), shape=process, colour=process)) +
  geom_point() + scale_shape_manual(values=c(5,6,17,19)) +
  geom_vline(xintercept=1, linetype="dashed") + theme_bw() +
  theme(legend.key.height=unit(4, "mm"), legend.title=element_blank(),
        legend.key=element_blank(), axis.title.y=element_blank(),
        #panel.border=element_rect(colour="black"), 
        panel.grid.major=element_blank()) +
  scale_x_log10(name="Retention factor")
  scale_x_continuous(name="Retention factor", trans=log_trans(), breaks=base_breaks(), labels=prettyNum) +
  scale_colour_brewer(palette="Dark2") + #ggtitle("A")
  
ggsave("black bean slide.png", width=180, height=100, unit="mm", dpi=300, plot=fig6a)

#with lattice
library(latticeExtra)
plot <- xyplot(reorder(compound, mean.rf) ~ mean.rf, data=bb, groups=process,
               pch=c(1,2,3,4), xlab="Mean retention factor", ylab=NULL, auto.key=T,
               scales=list(x=list(log=10, equispaced.log=FALSE)))

trellis.device(device="png", filename="black bean slide lattice.png")
print(plot)
dev.off()

#--------------------------------------------------------------------------------------

#RFs of quercetin, Q3Rut, 5CQA and caffeic acid, on log scales
#Quercetin
q <- read.csv("q_only2.csv", header=T)
#q <- rfs %>% filter(compound == "Quercetin")
fig5a <- ggplot(q, aes(x=mean.rf, y=reorder(food2, mean.rf), shape=process, colour=process)) + 
  geom_point() + geom_vline(xintercept=1, linetype="dashed") +
  theme(axis.title=element_blank(), panel.border=element_rect(colour="black"), 
        legend.title=element_blank(), legend.key.height=unit(4, "mm"),
        panel.grid.major=element_blank(), legend.key=element_blank()) + #ggtitle("A") +
  scale_colour_brewer(palette="Dark2") +
  scale_shape_manual(values=c(21,16,15,17,14,6,5,19)) + scale_x_log10(breaks=c(0.2,0.5,1,2))

ggsave("Quercetin for slide.png", width=140, height=70, unit="mm", dpi=600, plot=fig5a)

#Q3Rut
q3r <- read.csv("q3r_only2.csv")
fig5b <- ggplot(q3r, aes(x=mean.rf, y=reorder(food2, mean.rf), shape=process, colour=process)) +
  geom_point() + 
  theme(axis.title.y=element_blank(), panel.grid.major=element_blank(), legend.key.height=unit(4, "mm"),
        legend.title=element_blank(), #legend.position="none",
        panel.border=element_rect(colour="black"), legend.key=element_blank()) +
  scale_shape_manual(values=c(15,17,14,6,19)) + #ggtitle("B") +
  scale_x_log10(name="Retention factor", breaks=c(0.1,0.5,1,1.5)) +
  scale_colour_brewer(palette="Dark2") +
  geom_vline(xintercept=1, linetype="dashed")

ggsave("Q3Rut for slide.png", width=110, height=55, unit="mm", dpi=600, plot=fig5b)

#5-caffeoylquinic acid
fig4a <- ggplot(cqa, aes(x=mean.rf, y=reorder(food2, mean.rf), shape=process)) + 
  geom_point() + geom_vline(xintercept=1, linetype="dashed") +
  theme(legend.position="none", panel.border=element_rect(colour="black"),
        axis.title = element_blank(), panel.grid.major=element_blank(), legend.title=element_blank()) +
  scale_x_log10(breaks=c(0.05,0.1,0.2,0.5,1,2)) + ggtitle("A") +
  scale_shape_manual(values=c(21,15,17,14,18,6,5,13,19))

ggsave("fig4a revised.png", width=80, height=60, unit="mm", dpi=600, plot=fig4a)

#caffeic acid RFs on a log scale
ca <- read.csv("ca_only2.csv", header=T)
fig4b <- ggplot(ca, aes(x=mean.rf, y=reorder(food2, mean.rf), shape=process)) +
  geom_point() + geom_vline(xintercept=1, linetype="dashed") +
  theme(legend.key.height=unit(4, "mm"), panel.border=element_rect(colour="black"), 
        axis.title.y = element_blank(), legend.title=element_blank(),
        panel.grid.major=element_blank(), legend.key=element_blank()) +
  scale_x_log10(name="Retention factor", breaks=c(0.2,1,0.5,1,2,5)) +
  scale_shape_manual(values=c(21,16,15,17,14,18,6,5,13,19)) + ggtitle("B")

ggsave("fig4b revised.png", width=120, height=70, unit="mm", dpi=600, plot=fig4b)

#---------------------------------------------------------------------------------------

#potato scatter
potato <- read.csv("potato.csv", header=T)
#may need to find location of file to source
source("base_breaks.R")

fig6b <- ggplot(potato, aes(x=mean.rf, y=reorder(compound, mean.rf), shape=process)) +
  geom_point() + scale_shape_manual(values=c(21,5,15,18,19)) + 
  geom_vline(xintercept=1, linetype="dashed") +
  ggtitle("B") +
  theme(legend.key.height=unit(4, "mm"), axis.title.y=element_blank(), 
        panel.border=element_rect(colour="black"), legend.title=element_blank(),
        panel.grid.major=element_blank(), legend.key=element_blank()) +
  scale_x_continuous(name="Retention factor", trans=log_trans(),
                     breaks=base_breaks(), labels=prettyNum)

ggsave("fig6b revised.png", width=171, height=50, unit="mm", dpi=600, plot=fig6b)

#--------------------------------------------------------------------------------------------

#facet scatter RFs folin. subset Folin data only
rfs.folin <- subset(rfs, method == "Folin assay")
rfs.folin1 <- rfs %>% filter(method == "Folin assay", process == "Boiled")
ggplot(rfs.folin, aes(x=process, y=mean.rf)) + geom_point() + facet_grid(food.group ~ .) + theme_bw() +
  geom_jitter(position=position_jitter(w=0.2, h=0.15)) + scale_y_continuous(name="Mean retention factor") +
  scale_x_discrete(limits=c("Baked", "Blanched", "Boiled", "Fried", "Grilled",
                            "Microw.", "Pr. boiled", "Pr. steamed", "Steamed"))

#with lattice
library(latticeExtra)
xyplot(process ~ mean.rf, groups=food.group, data=rfs.folin)

#--------------------------------------------------------------------------------------------

#scatter RFs by process
ggplot(rfs, aes(x=process, y=mean.rf, colour=compound.class)) +
  geom_jitter(position=position_jitter(w=0.2, h=0.15)) + 
  #geom_hline(yintercept=1, linetype="dashed") + 
  #annotate("text", x=0.95, y=1.2, label="Mean RF=1") +
  theme_bw() + labs(x="Process", y="Aggregated retention factor") + 
  scale_x_discrete(limits=c("Baked", "Blanched", "Boiled", "Fried", "Grilled",
                            "Microw.", "Pr. boiled", "Pr. steamed", "Steamed")) +
  scale_y_log10(name="Mean retention factor") + theme_bw() +
  theme(legend.position=c(0.13,0.8))

ggsave("scatter rfs by process.png", width=6.5, height=4, unit="in", plot=processes)
                      
#--------------------------------------------------------------------------------------

# ANOVA RFs (experimental, not used in publication)

library(dplyr)
rfsfolin <- rfs %>% 
  filter(method=="Folin assay") %>% 
  select(food, process, mean.rf) %>%
  mutate(log.rf = log2(mean.rf))

summary(rfsfolin)

#plot observations
plot(rfsfolin$mean.rf)
plot(rfsfolin$log.rf)
qqnorm(rfsfolin$log.rf) # distribution not far from normal

#two way ANOVA with food and process
rfs.aov <- aov(log.rf ~ food * process, data=rfsfolin) 
summary(rfs.aov)

#diagnostic plots
plot(rfs.aov, which=c(1:2))
#too many factor levels -> model not valid???

############################################################################

rfschrom <- rfs %>% 
  filter(method=="Chrom.") %>% 
  select(food, process, compound, mean.rf) %>%
  mutate(log.rf = log2(mean.rf))

summary(rfschrom)

#-----------------------------------------------------------------------------------

# Barplots RF counts

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

#---------------------------------------------------------------------------------

# Boxplots RFs

library(ggplot2)

#Will's boxplot

g <- ggplot(boxplot.df, aes(x=Box, y=XObs_gLog_QC.LSC)) +
  geom_boxplot(outlier.colour="red", outlier.size=0) +
  stat_summary(fun.y=mean, geom="point", shape=18, size=12, col="red", fill="red") +
  geom_jitter(position=position_jitter(w=0.15,h=0.15)) +
  geom_hline(yintercept=b1$stats[[5]], colour="red", size=1) +
  annotate("text", label=Above.zero.samples, y=max(nonzero.boxplot[,j])*1.1, x=2, size=3.5) +
  annotate("text", label=Zero.outliers, y=max(nonzero.boxplot[,j])*1.1,x=1,size=3.5) +
  labs(title=Titlename) + theme_bw(20)

#Box and points 16 processes

rfs <- read.csv("raw_data.csv", header=T)

fig3 <- ggplot(rfs, aes(x=process.full, y=mean.rf)) +
  geom_boxplot(outlier.colour="red",outlier.size=0) +
  geom_jitter(position=position_jitter(w=0.15, h=0), size=2, colour="grey60", alpha=.6) +
  #stat_summary(fun.y=mean, geom="point",shape=18, size=4) +
  geom_hline(yintercept=1, linetype="dashed") + coord_flip() +
  theme_bw(base_size=9) + 
  theme(panel.border=element_rect(colour="black"), axis.title.y=element_blank(),
        panel.grid.major=element_blank()) +
  scale_y_continuous(name="Retention factor", trans=log_trans(), breaks=base_breaks(), labels=prettyNum) +
  scale_x_discrete(limits=c("Pressure-boiled", "Pressure-steamed", "Microwaved", "Canned",
                            "Steamed", "Boiled", "Fried", "Frozen", "Blanched", "Jam making",
                            "Stored at room temperature", "Frozen, stored frozen",
                            "Stored refrigerated", "Stored frozen", "Pasteurized", "Grilled"),
                   labels=c("Pressure\nboiled", "Pressure\nsteamed", "Microwaved", "Canned",
                            "Steamed", "Boiled", "Fried", "Frozen", "Blanched", "Jam-\nmaking",
                            "Stored at\nroom temp", "Frozen,\nstored frozen",
                            "Stored\nrefrigerated", "Stored\nfrozen", "Pasteurized", "Grilled"))

ggsave("fig3revised.png", width=80, height=110, units="mm", dpi=600, plot=fig3)

#in lattice
library(lattice)
bwplot(reorder(process, mean.rf, FUN=median) ~ mean.rf, data=rfs, jitter.data=T,
       xlab="Mean retention factor",
       scales=list(x=list(log=10, equispaced.log=F)), subset=rfs$top.process=="yes")

#----------------------------------------------------------------------------------------------------

#Boxplot all RFs by method
ggplot(rfs, aes(x=method, y=mean.rf)) + geom_boxplot(outlier.size=1.5, varwidth=T) + theme_bw() +
  geom_hline(yintercept=1, linetype="dashed") + theme(axis.title.x=element_blank()) +
  scale_y_log10(name="Retention factor", breaks=c(0.01,0.05,0.1,0.5,1,5,10,50)) + 
  scale_x_discrete(labels=c("Chromatography", "Chromatography\nwith\nhydrolysis",
                            "Folin-\nCiocalteu\nassay", "pH\ndifferential\nassay"))

ggsave("Boxplot all RFs by method.png", width=133, height=100, unit="mm")

#with base graphics and nice colours
library(RColorBrewer)
boxplot(log(rfs$mean.rf) ~ method, data = rfs, xlab="Retention factor", horizontal = T, col=brewer.pal(4,"Set3"))

#Frequency polygon
ggplot(rfs, aes(x=nraw.rf)) + geom_freqpoly(colour="black") + theme_bw() + 
  labs(x="Raw retention factors", y="Aggregated retention factors") 
ggsave("aggregation histogram.png", width=133, height=100, unit="mm")

#Boxplot RFs by method
ggplot(rfs, aes(x=method, y=mean.rf)) + geom_boxplot(outlier.size=1.5) + theme_bw() +
  geom_hline(yintercept=1, linetype="dashed") + labs(x="", y="Retention factor") +  
  scale_y_continuous(expand = c(0.02,0.02), limits=c(0,7.5), breaks=c(0:8))

ggsave("Boxplot all RFs by method.png", width=4, height=3, unit="in")
                      