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
