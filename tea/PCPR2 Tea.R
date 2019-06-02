# Load the data in the X matrix containing NMR spectra and Z matrix containing the list of explanatory variables of interest
library(pcpr2)

tea <- read.csv("data/Tea mod joe.csv")

# Subset X and Y variables
X_DataMatrix <- log2(tea[, 6:72]) %>% as.matrix
Y_var <- tea[, 2:4]
Y_var <- Y_var %>% mutate_at(vars(Replicate, Experiment), as.factor)

runPCPR2(X_DataMatrix, Y_var) %>% plotProp