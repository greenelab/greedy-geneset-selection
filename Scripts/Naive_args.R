# |r|, redundancy, and DM size
args <- commandArgs(trailingOnly=TRUE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LIBRARIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#load required libraries
library(ff)
library(doMC)
#library(multicore)
#library(fastcluster)
library(curatedOvarianData)

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
    ncore = detectCores()
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
  data.correlation <- bigcorPar(data, nblocks=nb)
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

getEligibleGenes <- function(corBin, rP, redundancy)
{
  
  values <- lapply(rP, function(x){
    corBinMat <- corrMat2BinMat(corBin, x)
    return(sum(colSums(corBinMat) > redundancy))
  })
  return(unlist(values))
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ANALYSIS OF THE GENE EXPRESSION DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#===================================
#CHECK ARGUMENTS
#===================================
if(args[1] > 1 | args[1] < 0)
{
  stop("The correlation threshold must be between 0 and 1")
  
}
if(args[2] < 1 | as.numeric(args[2]) %% 1 != 0)
{
  stop("Redundancy must be a positive integer")
}


#===================================
#CREATE BINARY MATRICES OF TCGA DATA
#===================================

#get the good sample IDs
tcga.samples <- read.csv("Data/SampleFilter/TCGA_eset_samplesRemoved_noDoppel.csv")


data(TCGA_eset)
TCGA <- exprs(TCGA_eset)
TCGA <- TCGA[,as.character(tcga.samples$x)]
rm(TCGA_eset)
gc()


TCGA.quantile <- apply(TCGA,1,function(x){
  return(quantile(x,0.90,na.rm=T))
})

TCGA <- TCGA[TCGA.quantile > 5,]

TCGA.cor <- getCorrelationMatrix(t(TCGA))

#create the binary adjacency matrix from the correlation matrix
TCGA.cor.bin <- corrMat2BinMat(TCGA.cor, as.numeric(args[1]))

#Get the number of edges per gene
edges <- colSums(TCGA.cor.bin)

#Sort genes by the number of edges
DM <- colnames(TCGA.cor.bin)[order(edges, decreasing=TRUE)]
DM <- DM[1:as.numeric(args[3])]

#Calculate the number of predictable genes
DM.matrix <- TCGA.cor.bin[,DM]

#Which rows have at least redundancy edges
predictable <- rownames(DM.matrix)[rowSums(DM.matrix) > as.numeric(args[2])]
predictable <- setdiff(predictable, DM)


print("Directly Measured Genes:")
print("-----------------------------------------------------------------------------------------------")
print("-----------------------------------------------------------------------------------------------")
cat("Directly Measured Genes: ")
print(paste(DM, collapse=", "))
print("-----------------------------------------------------------------------------------------------")
print("-----------------------------------------------------------------------------------------------")
cat("Predictable Genes: ")
print(paste(predictable, collapse=", "))
print("-----------------------------------------------------------------------------------------------")
print("-----------------------------------------------------------------------------------------------")
print(paste("Score:", length(predictable), sep=" "))


