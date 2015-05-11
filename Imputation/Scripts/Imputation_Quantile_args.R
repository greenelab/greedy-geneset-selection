args <- commandArgs(trailingOnly=TRUE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LIBRARIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#library(ggplot2)
#library(gplots)
library(ff)
library(doMC)
#library(multicore)
#library(fastcluster)
#library(hexbin)
library(curatedOvarianData)
library(randomForest)
library(survcomp)
#source("http://bioconductor.org/biocLite.R")
#biocLite("survcomp")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#CONSTANTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SEED <- 0
set.seed(SEED)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#FUNCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#The following function accelerates the creation of large correlation 
#matrices. Function taken from:
#https://gist.github.com/bobthecat/5024079
#http://www.r-bloggers.com/large-correlation-in-parallel/

bigcorPar <- function(x, nblocks = 10, verbose = TRUE, ncore="all", ...){
  library(ff, quietly = TRUE)
  require(doMC)
  if(ncore=="all"){
    ncore = multicore:::detectCores()
    registerDoMC(cores = ncore)
  } else{
    registerDoMC(cores = ncore)
  }
  
  cat("Using ", ncore, " cores\n")
  NCOL <- ncol(x)
  
  ## test if ncol(x) %% nblocks gives remainder 0
  if (NCOL %% nblocks != 0){stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")}
  
  ## preallocate square matrix of dimension
  ## ncol(x) in 'ff' single format
  corMAT <- ff(vmode = "single", dim = c(NCOL, NCOL))
  cat("Allocated Memory for the Resulting Matrix\n")
  
  ## split column numbers into 'nblocks' groups
  SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))
  cat("Split the Data into Blocks\n")
  
  
  ## create all unique combinations of blocks
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)
  cat("Calculate Combinations of Blocks")
  
  ## iterate through each block combination, calculate correlation matrix
  ## between blocks and store them in the preallocated matrix on both
  ## symmetric sides of the diagonal
  results <- foreach(i = 1:nrow(COMBS)) %dopar% {
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")
    flush.console()
    COR <- cor(x[, G1], x[, G2], use="complete.obs", ...)
    corMAT[G1, G2] <- COR
    corMAT[G2, G1] <- t(COR)
    COR <- NULL
  }
  
  gc()
  return(corMAT)
}

#Get the number of cpu blocks to use when creating the correlation matrix
#The total number of genes should be evenly divided by the number of blocks
getNumBlocks <- function(nGenes)
{
  rs <- 1
  
  for(i in 2:20)
  {
    if( nGenes %% i == 0)  
    {
      rs <- i
    }
  }
  return(rs)
}


#calculate a correlation for each column to each other column
getCorrelationMatrix <- function(data)
{
  nb <- getNumBlocks(ncol(data))
  data.correlation <- bigcorPar(data, nblocks=nb, ncore=4)
  correlation.matrix <- as.matrix(data.correlation[,])
  colnames(correlation.matrix) <- colnames(data)
  rownames(correlation.matrix) <- colnames(data)
  return(correlation.matrix)
}

#convert a continuous correlation matrix into a 0/1 representation based on greater than/equal to a threshold
corrMat2BinMat <- function(mat, thresh)
{
  
  result <- (abs(mat) >= thresh) * 1
  return(result)
}

MAD <- function(x)
{
  thisMedian <- median(x, na.rm=TRUE)
  result <- median(abs(x - thisMedian), na.rm=TRUE)
  return(result)
}

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
#ANALYSIS OF THE GENE EXPRESSION DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


data(TCGA_eset)
TCGA <- exprs(TCGA_eset)
rm(TCGA_eset)
gc()
#load the good samples
tcga.samples <- read.csv("Data/SampleFilter/TCGA_eset_samplesRemoved_noDoppel.csv")
TCGA <- TCGA[,as.character(tcga.samples$x)]
#permute the samples
TCGA <- TCGA[,sample(ncol(TCGA))]
sampleTrainingUpperBound <- round(ncol(TCGA) * (1/3) * 2)
sampleTestingLowerBound <- sampleTrainingUpperBound + 1
sampleTestingUpperBound <- ncol(TCGA)

#Let's test in Tothill Data
data(GSE9891_eset) 
Tothill <- exprs(GSE9891_eset)
rm(GSE9891_eset)
gc()
#load the good samples
tothill.samples <- read.csv("Data/SampleFilter/GSE9891_eset_samplesRemoved_noDoppel.csv")
Tothill <- Tothill[,as.character(tothill.samples$x)]

#Let's test in Yoshihara Data
data(GSE32062.GPL6480_eset)
Yoshihara <- exprs(GSE32062.GPL6480_eset)
rm(GSE32062.GPL6480_eset)
gc()
#load the good samples
yoshihara.samples <- read.csv("Data/SampleFilter/GSE32062.GPL6480_eset_samplesRemoved_noDoppel.csv")
Yoshihara <- Yoshihara[,as.character(yoshihara.samples$x)]

#Let's test in the Bonome Data (185 Samples from USA; Memorial Sloan-Kettering)
data(GSE26712_eset)
#remove the 10 normal samples from the Bonome data
Bonome <- exprs(GSE26712_eset)
bonome.samples <- read.csv("Data/SampleFilter/GSE26712_eset_samplesRemoved_noDoppel.csv")
Bonome <- Bonome[,as.character(bonome.samples$x)]
#Bonome <- Bonome[,GSE26712_eset$sample_type == "tumor"]
rm(GSE26712_eset)
gc()

#Let's get the Crijns data
data(GSE13876_eset)
Crijns <- exprs(GSE13876_eset)
Crijns.samples <- read.csv("Data/SampleFilter/GSE13876_eset_samplesRemoved_noDoppel.csv")
Crijns <- Crijns[,as.character(Crijns.samples$x)]
rm(GSE13876_eset, Crijns.samples)
gc()

# #Let's get the Dressman data
# data(PMID17290060_eset)
# Dressman <- exprs(PMID17290060_eset)
# Dressman.samples <- read.csv("Data/SampleFilter/PMID17290060_eset_samplesRemoved.csv")
# Dressman <- Dressman[,as.character(Dressman.samples$x)]
# rm(PMID17290060_eset, Dressman.samples)
# 
# #Let's get the Pils data
# data(GSE49997_eset)
# Pils <- exprs(GSE49997_eset)
# Pils.samples <- read.csv("Data/SampleFilter/GSE49997_eset_samplesRemoved.csv")
# Pils <- Pils[,as.character(Pils.samples$x)]
# rm(GSE49997_eset, Pils.samples)
# gc()
# 
# 
# #Let's get the Karlan data
# data(GSE51088_eset)
# Karlan <- exprs(GSE51088_eset)
# Karlan.samples <- read.csv("Data/SampleFilter/GSE51088_eset_samplesRemoved.csv")
# Karlan <- Karlan[,as.character(Karlan.samples$x)]
# rm(GSE51088_eset, Karlan.samples)
# gc()


#Let's get the TCGA RNAseq data
data(TCGA.RNASeqV2_eset)
TCGA.RNA <- exprs(TCGA.RNASeqV2_eset)
tcga.rna.samples <- intersect(tcga.samples$x, colnames(TCGA.RNA))
TCGA.RNA <- TCGA.RNA[,as.character(tcga.rna.samples)]
# common.genes <- intersect(rownames(TCGA), rownames(TCGA.RNA))
# common.genes.cor <- lapply(common.genes, function(x){
#   return(cor(TCGA[x,colnames(TCGA.RNA)], TCGA.RNA[x,]))
# })
# hist(unlist(common.genes.cor), main="Correlation of RNAseq and Affymetrix gene expression", xlab="Pearson's Correlation Coefficient")
# 
# TCGA.quantile <- unlist(apply(TCGA, 1, function(x){
#   quantile(x, prob=0.90)
# }))
# 
# genes.expressed <- rownames(TCGA)[TCGA.quantile > 5]
# common.genes.expressed <- intersect(genes.expressed, rownames(TCGA.RNA))
# common.genes.expressed.cor <- lapply(common.genes.expressed, function(x){
#   return(cor(TCGA[x,colnames(TCGA.RNA)], TCGA.RNA[x,]))
# })
# hist(unlist(common.genes.expressed.cor), main="Correlation of Quantile filtered\nRNAseq and Affymetrix gene expression", xlab="Pearson's Correlation Coefficient")
rm(TCGA.RNASeqV2_eset)
gc()

#Let's get the Goode data (384 samples measured on an Agilent wholegenome_4x44k_v1 platform
Goode <- read.table("Data/jen.comb.mtx2_exprs_Normalizer.txt", sep="\t")
Goode <- Goode[,-1]
mayo.samples <- read.csv("Data/SampleFilter/mayo.eset_samplesRemoved_noDoppel.csv")
Goode <- Goode[,as.character(mayo.samples$x)]
cat("Loaded Goode Data \n")



print("Generating Correlation Matrix...")
TCGA.cor <- getCorrelationMatrix(t(TCGA))

print(paste("Running the imputation experiment for each parameter combination...(", length(args), ")", sep=""))
for(k in 1:length(args))
{

  print(args[k])
  corThresh <- as.numeric(paste("0.",tail(unlist(strsplit(args[k], "\\.")), n=1), sep=""))
  TCGA.cor.bin <- corrMat2BinMat(TCGA.cor, corThresh)

  measured.genes <- c()
  covered.genes <- c()
  prefix <- c()

  args.detailed <- unlist(strsplit(args[k], "\\."))
  if(length(args.detailed) == 6)
  {
    candidateName <- args.detailed[1]
    #print(paste("This uses the ", candidateName, " candidate gene list", sep=""))
    prefix <- c(prefix,paste(candidateName, ".Candidate.Genes", sep=""))
    thisDir <- paste("Data/Quantile.GGS.Parameter.Sweep.Results/", candidateName, ".Candidate.Genes/genesets/", sep="")

    fname.cleaned <- paste(args.detailed[2:6], collapse=".")

    measured.fname <- paste(thisDir,fname.cleaned,".measured.txt",sep="")
    measured.genes <- read.csv2(measured.fname, header=F, sep=",")
    measured.genes <- as.vector(unlist(measured.genes))
    #setdiff(measured.genes, colnames(TCGA.cor.bin))

    #print(measured.fname)

    covered.fname <- paste(thisDir,fname.cleaned,".covered.txt",sep="")
    covered.genes <- read.csv2(covered.fname, header=F, sep=",")
    covered.genes <- as.vector(unlist(covered.genes[1,]))
    covered.genes <- setdiff(covered.genes, measured.genes)
  }
  
  else
  {
    measured.fname <- paste("Data/Quantile.GGS.Parameter.Sweep.Results/genesets/",args[k],".measured.txt",sep="")
    measured.genes <- read.csv2(measured.fname, header=F, sep=",")
    measured.genes <- as.vector(unlist(measured.genes))
    setdiff(measured.genes, colnames(TCGA.cor.bin))


    covered.fname <- paste("Data/Quantile.GGS.Parameter.Sweep.Results/genesets/",args[k],".covered.txt",sep="")
    covered.genes <- read.csv2(covered.fname, header=F, sep=",")
    covered.genes <- as.vector(unlist(covered.genes[1,]))
    covered.genes <- setdiff(covered.genes, measured.genes)
  }
  
  print(setdiff(measured.genes, rownames(TCGA)))
  TCGA.measured <- normalize.JR(TCGA[measured.genes,])
  TCGA.estimated <-normalize.JR(TCGA[covered.genes,])

  print(paste("#Measured Genes: ", length(measured.genes), "\t#Covered Genes: ", length(covered.genes)))
  #Subset to only the measured genes
  Tothill.measured <- normalize.JR(Tothill[measured.genes,])
  Yoshihara.measured <- normalize.JR(Yoshihara[intersect(measured.genes, rownames(Yoshihara)),])
  Bonome.measured <- normalize.JR(Bonome[measured.genes,])
  Crijns.measured <- normalize.JR(Crijns[intersect(measured.genes, rownames(Crijns)),])
#   Dressman.measured <- normalize.JR(Dressman[intersect(measured.genes, rownames(Dressman)),])
#   Pils.measured <- normalize.JR(Pils[intersect(measured.genes, rownames(Pils)),])
#   Karlan.measured <- normalize.JR(Karlan[intersect(measured.genes, rownames(Karlan)),])
  Goode.measured <- normalize.JR(Goode[intersect(measured.genes, rownames(Goode)),])
  TCGA.RNA.measured <- normalize.JR(TCGA.RNA[intersect(measured.genes, rownames(TCGA.RNA)),])

  
  #The empty object which will eventually hold our result data frame
  result <- c()
  
  #Make sure that there really are covered genes
  if(length(covered.genes[!is.na(covered.genes)]) < 1)
  {
    fname <- paste("Imputation/Data/",args[k],".Quantile.imputation.results.txt", sep="")
    write.table(c(), fname, quote=FALSE)
    next
  }
  #Loop over the genes in the covered set. For each gene we need to perform an imputation experiment
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
      model <- randomForest(x=TCGA.measured[1:sampleTrainingUpperBound,tags], y=TCGA.estimated[1:sampleTrainingUpperBound,thisGene])
      
      #Test the model in the remaining 1/3 of the TCGA Data
      #We perform a prediction and correlate it to the ground truth expression.
      prediction.TCGA <- predict(model, TCGA.measured[sampleTestingLowerBound:sampleTestingUpperBound,tags], type="response")
      groundTruth.TCGA <- normalize.Vec(TCGA.estimated[sampleTestingLowerBound:sampleTestingUpperBound, thisGene])
      prediction.TCGA.cor.pearson <- cor(prediction.TCGA, groundTruth.TCGA)
      prediction.TCGA <- prediction.TCGA.cor.spearman <- cor(prediction.TCGA, groundTruth.TCGA, method="spearman")
    }
    else
    {
      #If we only have one predictor variable, then we use polynomial regression
      y <- TCGA.estimated[1:sampleTrainingUpperBound,thisGene]
      x <- TCGA.measured[1:sampleTrainingUpperBound,tags]
      model <- lm(y ~ poly(x, 3, raw=TRUE))
      
      
      prediction.TCGA <- predict(model, data.frame(x=TCGA.measured[sampleTestingLowerBound:sampleTestingUpperBound,tags]))
      groundTruth.TCGA <- normalize.Vec(TCGA.estimated[sampleTestingLowerBound:sampleTestingUpperBound, thisGene])
      prediction.TCGA.cor.pearson <- cor(prediction.TCGA, groundTruth.TCGA)
      prediction.TCGA <- prediction.TCGA.cor.spearman <- cor(prediction.TCGA, groundTruth.TCGA, method="spearman")
    }
    
    #Test in the Tothill Data
    prediction.Tothill <- predictExpression(tags, thisGene, Tothill, model, Tothill.measured)
    
    #Test in the Bonome Data
    prediction.Bonome <- predictExpression(tags, thisGene, Bonome, model, Bonome.measured)
    
    #We need to test whether the covered and 'tag' genes are present in the Goode dataset
    prediction.Goode <- predictExpression(tags, thisGene, Goode, model, Goode.measured)
    
    #We need to test whether the covered and 'tag' genes are present in the Yoshihara dataset
    prediction.Yoshihara <- predictExpression(tags, thisGene, Yoshihara, model, Yoshihara.measured)
    
    #We need to test whether the covered and 'tag' genes are present in the Crijns dataset
    prediction.Crijns <- predictExpression(tags, thisGene, Crijns, model, Crijns.measured)
    
    #We need to test whether the covered and 'tag' genes are present in the TCGA RNAseq dataset
    prediction.TCGA.RNA <- predictExpression(tags, thisGene, TCGA.RNA, model, TCGA.RNA.measured)
    
#     #We need to test whether the covered and 'tag' genes are present in the Dressman dataset
#     prediction.Dressman <- predictExpression(tags, thisGene, Dressman, model, Dressman.measured)
#     
#     #We need to test whether the covered and 'tag' genes are present in the Pils dataset
#     prediction.Pils <- predictExpression(tags, thisGene, Pils, model, Pils.measured)
#     
#     #We need to test whether the covered and 'tag' genes are present in the Pils dataset
#     prediction.Karlan <- predictExpression(tags, thisGene, Karlan, model, Karlan.measured)
    
    result <- rbind(result, c(thisGene, length(tags), prediction.TCGA, prediction.Tothill, prediction.Yoshihara, prediction.Bonome, prediction.Crijns, prediction.Goode, prediction.TCGA.RNA))
    progress <- paste((i / length(covered.genes) ) * 100, "%    \r")
    cat(progress)
  }
  
  result <- data.frame(result)
  rownames(result) <- result[,1]
  result <- result[,-1]
  colnames(result) <- c("nTagGenes", "TCGA.Spearman", "Tothill.Spearman", "Yoshihara.Spearman",
                        "Bonome.Spearman","Crijns.Spearman", "Goode.Spearman", "TCGA.RNAseq.Spearman")
  
#   result.mod <- data.frame(CorrelationCoefficients=result$TCGA.Spearman, CorrType="Spearman", DataSet="TCGA")
#   result.mod <- rbind(result.mod, data.frame(CorrelationCoefficients=result$Tothill.Spearman, CorrType="Spearman", DataSet="Tothill"))
#   result.mod <- rbind(result.mod, data.frame(CorrelationCoefficients=result$Yoshihara.Spearman, CorrType="Spearman", DataSet="Yoshihara"))
#   result.mod <- rbind(result.mod, data.frame(CorrelationCoefficients=result$Bonome.Spearman, CorrType="Spearman", DataSet="Bonome"))
#   result.mod <- rbind(result.mod, data.frame(CorrelationCoefficients=result$Crijns.Spearman, CorrType="Spearman", DataSet="Crijns"))
#   result.mod <- rbind(result.mod, data.frame(CorrelationCoefficients=result$Dressman.Spearman, CorrType="Spearman", DataSet="Dressman"))
#   result.mod <- rbind(result.mod, data.frame(CorrelationCoefficients=result$Goode.Spearman, CorrType="Spearman", DataSet="Goode"))

  
 

#   #Let's make sure that no dataset is empty
#   tbl <- table(result.mod$DataSet, !is.na(result.mod$CorrelationCoefficients))
#   tbl.empty <- rownames(tbl)[tbl[,2] == 0]
#   tbl.full <- setdiff(rownames(tbl), tbl.empty)
#   result.mod <- result.mod[result.mod$DataSet %in% tbl.full,]
#   
#   result.mod$CorrelationCoefficients <- as.numeric(paste(result.mod$CorrelationCoefficients))
#   fname <- paste("Imputation/Figures/",args[k],".Quantile.imputation.correlation.accuracy.png", sep="")
#   
#   ggplot(result.mod, aes(CorrType, CorrelationCoefficients)) + geom_boxplot() + facet_wrap(~DataSet)
#   ggsave(fname, width=11, height=8)
  
  
  #Save the results dataframe 
  fname <- paste("Imputation/Data/",args[k],".Quantile.imputation.results.txt", sep="")
  write.table(result, fname, quote=FALSE)

  if(FALSE)
  {
    if(tbl["TCGA",2] > 0)
    {
      png(paste("Imputation/Figures/",args[k],".MAD.imputation.correlation.accuracy.TCGA.testing.MAD.png", sep=""), width=1000, height=700)
      plot(result.var$MAD.TCGA, as.numeric(paste(result.var$TCGA.Spearman)), xlab="MAD", ylab="Imputation Spearman Rho", main="In TCGA Testing Partition\nAccuracy of Imputation by MAD")
      abline(v=0.3, lty=2, col="red")
      abline(v=0.5, lty=2, col="red")
      dev.off()
      
      png(paste("Imputation/Figures/",args[k],".MAD.imputation.correlation.accuracy.TCGA.testing.normally.expressed.MAD.png", sep=""), width=1000, height=700)
      plot(result.var$MAD.TCGA, as.numeric(paste(result.var$TCGA.Spearman)), col=(result.var$noExpress.TCGA>5)+1, xlab="MAD", ylab="Imputation Spearman Rho", main="In TCGA Testing Partition\nColored normally expressed")
      abline(v=0.3, lty=2, col="red")
      abline(v=0.5, lty=2, col="red")
      dev.off()
      
      png(paste("Imputation/Figures/",args[k],".MAD.imputation.correlation.accuracy.TCGA.testing.Expression.Level.png", sep=""), width=1000, height=700)
      plot(result.var$noExpress, as.numeric(paste(result.var$TCGA.Spearman)), col=(result.var$noExpress.TCGA>5)+1, xlab="99th Quantile Threshold", ylab="Imputation Spearman Rho", main="In TCGA Testing Partition\nColored normally expressed")
      dev.off()
      #plot(result.var$noExpress.TCGA, as.numeric(paste(result.var$TCGA.Spearman)), xlab="MAD", ylab="Imputation Spearman Rho", main="In TCGA Testing Partition\nAccuracy of Imputation by MAD")
      
    }
    
    if(tbl["Tothill",2] > 0)
    {
      png(paste("Imputation/Figures/",args[k],".MAD.imputation.correlation.accuracy.Tothill.MAD.png", sep=""), width=1000, height=700)
      plot(result.var$MAD.Tothill, as.numeric(paste(result.var$Tothill.Spearman)), xlab="MAD", ylab="Imputation Spearman Rho", main="In Tothill Data\nAccuracy of Imputation by MAD")
      dev.off()
    }
    
    if(tbl["Yoshihara",2] > 0)
    {
      png(paste("Imputation/Figures/",args[k],".MAD.imputation.correlation.accuracy.Yoshihara.MAD.png", sep=""), width=1000, height=700)
      plot(result.var$MAD.Yoshihara, as.numeric(paste(result.var$Yoshihara.Spearman)), xlab="MAD", ylab="Imputation Spearman Rho", main="In Yoshihara Data\nAccuracy of Imputation by MAD")
      dev.off()
    }
    
    if(tbl["Bonome",2] > 0)
    {
      png(paste("Imputation/Figures/",args[k],".MAD.imputation.correlation.accuracy.Bonome.MAD.png", sep=""), width=1000, height=700)
      plot(result.var$MAD.Bonome, as.numeric(paste(result.var$Bonome.Spearman)), xlab="MAD", ylab="Imputation Spearman Rho", main="In Bonome Data\nAccuracy of Imputation by MAD")
      dev.off()
    }
    
    if(tbl["Bentink",2] > 0)
    {
      png(paste("Imputation/Figures/",args[k],".MAD.imputation.correlation.accuracy.Bentink.MAD.png", sep=""), width=1000, height=700)
      plot(result.var$MAD.Bentink, as.numeric(paste(result.var$Bentink.Spearman)), xlab="MAD", ylab="Imputation Spearman Rho", main="In Bentink Data\nAccuracy of Imputation by MAD")
      dev.off()
      
      png(paste("Imputation/Figures/",args[k],".MAD.imputation.correlation.accuracy.Bentink.ReTransform.MAD.png", sep=""), width=1000, height=700)
      plot(result.var$MAD.Bentink, as.numeric(paste(result.var$Bentink.ReTransform.Spearman)), xlab="MAD", ylab="Imputation Spearman Rho", main="In Bentink Data\nAccuracy of Imputation by MAD")
      dev.off()
      
    }
    
    if(tbl["Goode",2] > 0)
    {
      png(paste("Imputation/Figures/",args[k],".MAD.imputation.correlation.accuracy.Goode.MAD.png", sep=""), width=1000, height=700)
      plot(result.var$MAD.Goode, as.numeric(paste(result.var$Goode.Spearman)), xlab="MAD", ylab="Imputation Spearman Rho", main="In Goode Data\nAccuracy of Imputation by MAD")
      dev.off()
    }
    
  }
 
}
























