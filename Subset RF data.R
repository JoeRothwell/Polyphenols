#to subset a data frame
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
                      
                      
                      