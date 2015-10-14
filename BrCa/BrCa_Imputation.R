#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LIBRARIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library("RTCGAToolbox")
library("randomForest")
library(curatedOvarianData)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#CONSTANTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GGS.results.dir <- "../Data/GGS.Parameter.Sweep.Results/genesets/"
correlation.thresholds <- c("0.60", "0.65", "0.70")
set.seed(0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#FUNCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

normalize.JR <- function(mydata)
{
  #taken from:
  #http://stackoverflow.com/questions/20046257/normalize-rows-of-a-matrix-within-range-0-and-1
  return(apply(mydata, 1, function(x){return(normalize.Vec(x))}))
}

normalize.Vec <- function(x)
{
  return((x-min(x, na.rm=TRUE))/(max(x, na.rm=TRUE)-min(x, na.rm=TRUE)))
}

predictExpression <- function(tags, thisGene, thisData, model, measured)
{
  tmp <- intersect(tags, rownames(thisData))
  if(length(tmp) == length(tags) & thisGene %in% rownames(thisData))
  {
    prediction <- NA
    measured <- measured[,tags]
    useable <- complete.cases(measured)
    if(length(tags) > 1)
    {
      
      prediction <- predict(model, measured[useable,], type="response") 
    }
    else
    {
      
      prediction <- predict(model,  data.frame(x=measured[useable]))
    }
    
    truth <- unlist(normalize.Vec(thisData[thisGene,useable]))
    prediction.spearman <- cor(prediction, truth, method="spearman")
    return(prediction.spearman)
  }
  else
  {
    return(NA)
  }
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#IMPUTATION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#Load the covered and predicted gene sets from the subfolders.
results.covered.file.list <- list.files(GGS.results.dir, pattern="*covered.txt", full.names=TRUE)
results.measured.file.list <- list.files(GGS.results.dir, pattern="*measured.txt", full.names=TRUE)

#Load the TCGA RNAseq breast cancer training/testing data set
BrCa.TCGA <- getFirehoseData("BRCA", runDate="20150402", RNAseq2_Gene_Norm = TRUE)
BrCa.RNAseq <- BrCa.TCGA@RNASeq2GeneNorm
rm(BrCa.TCGA)
gc()

#Randomly permute the samples
BrCa.RNAseq <- BrCa.RNAseq[,sample(ncol(BrCa.RNAseq), ncol(BrCa.RNAseq))]

#Divide samples into training (2/3) and testing (1/3) partitions
lastTrainingSample <- round(ncol(BrCa.RNAseq) * (2/3))
trainingID <- seq(from=1,to=lastTrainingSample)
testingID <- seq(from=lastTrainingSample+1,to=ncol(BrCa.RNAseq))

#Iterate over the correlation thresholds
for(corThresh in correlation.thresholds)
{
  
  TCGA.cor.bin <- read.delim(paste("BrCa.RNAseq.cor.", corThresh, ".txt", sep=""))
  
  subset.measured.file.list <- grep(corThresh, results.measured.file.list, value=TRUE) 
  subset.covered.file.list <- grep(corThresh, results.covered.file.list, value=TRUE)
  
  #Iterate over the paired measured and covered gene list files
  for(genelist.file.ID in 1:length(subset.measured.file.list))
  {
    
    
    #Initialize the covered and measured gene sets
    measured.fname <- subset.measured.file.list[genelist.file.ID]
    measured.genes <- read.csv2(measured.fname, header=F, sep=",", stringsAsFactors = FALSE)
    measured.genes <- paste(measured.genes)
    measured.genes <- intersect(measured.genes, rownames(BrCa.RNAseq))
    
    covered.fname <- subset.covered.file.list[genelist.file.ID]
    covered.genes <- read.csv2(covered.fname, header=F, sep=",", stringsAsFactors = FALSE)
    covered.genes <- paste(covered.genes)
    covered.genes <- setdiff(covered.genes, measured.genes)
    covered.genes <- intersect(covered.genes, rownames(BrCa.RNAseq))
    
    #Normalize the gene expression values for genes in the measured and covered sets
    TCGA.measured <- normalize.JR(BrCa.RNAseq[measured.genes,])
    TCGA.estimated <-normalize.JR(BrCa.RNAseq[covered.genes,])
    
    #Print status messages
    print(measured.fname)
    print(covered.fname)
    print(paste("#Measured Genes: ", length(measured.genes), "\t#Covered Genes: ", length(covered.genes)))
    
    #Initialize the result object
    result <- c()
    
    #For each predictable gene
    for(i in 1:length(covered.genes))
    {
      
      #Get the covered gene name
      thisGene <- covered.genes[i]
      
      #Identify which genes in the measured gene set 'tag' this covered gene
      neighbors <- rownames(TCGA.cor.bin)[TCGA.cor.bin[thisGene,] > 0]
      tags <- intersect(neighbors,measured.genes)
      
      model <- NA
      prediction.TCGA <- NA
      #This ensures that we use random forest regression only if there are 2 or more predictor variables.
      if(length(tags) > 1)
      {
        #Build the random forest regression model for this covered gene only using the measured genes that 'tag' it
        #The model is built only using 2/3 of the TCGA data.
        model <- randomForest(x=TCGA.measured[trainingID,tags], y=TCGA.estimated[trainingID,thisGene])
        
        #Test the model in the remaining 1/3 of the TCGA Data
        #We perform a prediction and correlate it to the ground truth expression.
        prediction.TCGA <- predict(model, TCGA.measured[testingID,tags], type="response")
        groundTruth.TCGA <- normalize.Vec(TCGA.estimated[testingID, thisGene])
        prediction.TCGA.cor.pearson <- cor(prediction.TCGA, groundTruth.TCGA)
        prediction.TCGA <- prediction.TCGA.cor.spearman <- cor(prediction.TCGA, groundTruth.TCGA, method="spearman")
      }
      
      else
      {
        #If we only have one predictor variable, then we use polynomial regression
        y <- TCGA.estimated[trainingID,thisGene]
        x <- TCGA.measured[trainingID,tags]
        model <- lm(y ~ poly(x, 3, raw=TRUE))
        
        
        prediction.TCGA <- predict(model, data.frame(x=TCGA.measured[testingID,tags]))
        groundTruth.TCGA <- normalize.Vec(TCGA.estimated[testingID, thisGene])
        prediction.TCGA.cor.pearson <- cor(prediction.TCGA, groundTruth.TCGA)
        prediction.TCGA <- prediction.TCGA.cor.spearman <- cor(prediction.TCGA, groundTruth.TCGA, method="spearman")
      }
      
      result <- rbind(result, c(thisGene, length(tags), prediction.TCGA))
      progress <- paste((i / length(covered.genes) ) * 100, "%    \r")
      cat(progress)
      
    }
    
    result <- data.frame(result)
    rownames(result) <- result[,1]
    result <- result[,-1]
    colnames(result) <- c("nTagGenes", "TCGA.Spearman")
    
    #Save the result data frame
    fname <- gsub(".covered.txt", ".imputation.results.txt", covered.fname)
    write.table(result, fname, quote=FALSE)
  }
  
  
}








