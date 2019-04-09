# Load the data in the X matrix containing NMR spectra and Z matrix containing the list of explanatory variables of interest
library(MetabolAnalyze)
library(car)

tea <- read.csv("Tea mod joe.csv")
X_DataMatrix <- log2(tea[, 6:72]) #X matrix
Z_InterestFactors <- tea[, 2:4] #Y variables

#myPath <- "Documents/PCPR2/" Metabolomics_data <- "X_MetaboMatrix.TXT"
#InterestFactors_data <- "Z_FactorMatrix.TXT"
#Metabo_FilePath = paste(myPath,Metabolomics_data, sep="")
#Factors_FilePath = paste(myPath,InterestFactors_data, sep="")
#X_DataMatrix <- read.delim(Metabo_FilePath, row.names = 1, header = TRUE, sep = "\t")
X_DataMatrixCentered = scale(X_DataMatrix, center = TRUE, scale = FALSE) 
X_DataMatrixScaled = scaling(X_DataMatrixCentered, type = "pareto")

# Center the data / Scale the data, edit the parameter "pareto" or "unit" of scaling according to your need
#Z_InterestFactors <- read.delim(Factors_FilePath, sep = "\t", header = TRUE, row.names = 1)
Z_InterestFactorsRowN <- nrow(Z_InterestFactors)
Z_InterestFactorsColN <- ncol(Z_InterestFactors) 
ColNames <- names(Z_InterestFactors)

# Obtain eigenvectors
pct_threshold = .95 # Amount of variability desired to be explained, to be edited with your preferences
X_DataMatrixScaled_transposed = t(X_DataMatrixScaled)
Mat2 <- X_DataMatrixScaled%*%X_DataMatrixScaled_transposed 
eigenData <- eigen(Mat2)
eigenValues = eigenData$values
ev_n <- length(eigenValues)
eigenVectorsMatrix = eigenData$vectors
eigenValuesSum = sum(eigenValues)
percents_PCs = eigenValues /eigenValuesSum
my_counter_2 = 0
my_sum_2 = 1
for (i in ev_n:1){
  my_sum_2 = my_sum_2 - percents_PCs[i]
  if ((my_sum_2) <= pct_threshold ){ 
    my_counter_2 = my_counter_2 + 1
  } 
}
if (my_counter_2 < 3){ pc_n =3
}else {
  pc_n = my_counter_2
}
pc_data_matrix <- matrix(data = 0, nrow = (Z_InterestFactorsRowN*pc_n), ncol = 1) 
mycounter = 0
for (i in 1:pc_n){
  for (j in 1:Z_InterestFactorsRowN){ 
    mycounter <- mycounter + 1                          
    pc_data_matrix[mycounter,1] = eigenVectorsMatrix[j,i] 
  }
}
AAA <- Z_InterestFactors[rep(1:Z_InterestFactorsRowN,pc_n),]
Data <- cbind(pc_data_matrix, AAA)
                                      
#Perform linear multiple regression models on each eigenvector with factors of interest as explanatory variables
#Categorical variables should be processed by as.factor, whereas continuous variables should not. 
#To be edited with your factors names

Z_InterestFactors$Replicate <- as.factor(Z_InterestFactors$Replicate)
Z_InterestFactors$Experiment <- as.factor(Z_InterestFactors$Experiment)
Z_InterestFactors$Altitude <- Z_InterestFactors$Altitude
DataCol <- ncol(Data)
typeIIIMatrix <- matrix(data = 0, nrow = (pc_n), ncol = (DataCol) ) 
ST_ResidualR2 <- matrix(data = 0, nrow = (pc_n), ncol = 2)                                      
for (i in 1:pc_n){                                        
  y = (((i-1)*Z_InterestFactorsRowN)+1)                                                                                                                      
  TotSumSq <- var(Data[y:(((i- 1)*Z_InterestFactorsRowN)+Z_InterestFactorsRowN),1])*(Z_InterestFactorsRowN-1)
                                                                               
#Edit the linear model with your factors
                                        
Model <- lm((pc_data_matrix ~ Replicate + Experiment + Altitude), Data[y:(((i-1)*Z_InterestFactorsRowN)+Z_InterestFactorsRowN),])
        AnalysisVariance <- Anova(Model,type=c(3))
        SumSq <- AnalysisVariance[1] 
        Residuals <- SumSq[DataCol+1,] 
        RR <- Residuals/TotSumSq
        R2 = 1 - RR
        ST_ResidualR2[i,] <- c(R2,RR)
        ST_ResidualR2_Names <- c("ST_R2", "ST_Residuals")
        colnames(ST_ResidualR2) = ST_ResidualR2_Names
        for (j in 1:(DataCol)){
        typeIIIMatrix[i,j] = as.numeric(SumSq[j+1,1])
        typeIIIMatrixNames <- c(ColNames, "SumSqResiduals")
        colnames(typeIIIMatrix) = typeIIIMatrixNames
    } 
}
partialR2Matrix <- matrix(data = 0, nrow = (pc_n), ncol = (DataCol-1) )
for (i in 1:pc_n){
    for (j in 1:(DataCol-1)){
      partialR2Matrix[i,j] = typeIIIMatrix[i,j] / (typeIIIMatrix[i,(DataCol)] + typeIIIMatrix[i,j]) 
    }
}
  partialR2MatrixWtProp <- matrix(data = 0, nrow = (pc_n), ncol = (DataCol)) 
for (i in 1:pc_n){
    weight = eigenValues[i]/sum(eigenValues[1:pc_n]) 
    for (j in 1:DataCol-1){
      partialR2MatrixWtProp[i,j] = partialR2Matrix[i,j]*weight
      partialR2MatrixWtProp[i,DataCol] = ST_ResidualR2[i,1]*weight 
    }
}
pR2Sums <- colSums(partialR2MatrixWtProp)*100 
plotnames = c( ColNames, "R2")
bp <- barplot(pR2Sums, xlab = "Factors", ylab = "Weighted Rpartial2", ylim= c(0,60),col = c("red"), las=2)
axis(1, at = bp, labels = plotnames, xlab = "Factors", cex.axis = 0.5, las=2) 
values = pR2Sums
new_values = round(values, 3)
text(bp, pR2Sums, labels = new_values, pos=3, cex = 0.8)

#output <- data.frame(plotnames, pR2Sums)
#library(ggplot2)
#ggplot(data=output, aes(x=plotnames, y=pR2Sums)) + geom_bar(stat="identity", colour="black", fill="red") +
  #theme_bw() + geom_text(aes(label=new_values),  vjust=-0.5) + ylim(0,60) +
  #scale_x_discrete(name="Factor", limits=c("Replicate", "Experiment", "Altitude", "R2"))

#library(lattice)
#barchart(pR2Sums~plotnames, data=output, origin=0, col="pink")
