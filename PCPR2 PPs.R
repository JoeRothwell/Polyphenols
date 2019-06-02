library(tidyverse)
library(zoo)

#read in X matrix and Y explanatory variables. Impute unknowns with medians
profilespp <- read.csv("data/Urinary polyphenols_EPIC.csv", header=T)

# Subset compound concentrations and metadata
mat <- profilespp %>% select(starts_with("uPPc"))
X_DataMatrix <- na.aggregate(mat, FUN = function(x) min(x)/2) %>% scale

# Y-variables. Remove country as directly correlated with centre
Z_InterestFactors <- profilespp %>% select(Batch:RE_ENERGY, -Country)

library(pcpr2)
runPCPR2(X_DataMatrix, Z_InterestFactors) %>% plotProp
