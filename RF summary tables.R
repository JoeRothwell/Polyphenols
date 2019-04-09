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

#------------------------------------------------------------------------------------------
#to summarise RF data by process with count of raw and aggregated RFs and median and IQR of aggregated RFs
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