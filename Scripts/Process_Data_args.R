
args <- commandArgs(trailingOnly=TRUE)
redundancy <- c(1,2,3)
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

tcga.candidates <- read.delim("Data/Candidate.Genes/TCGA2011.sig.mapped.csv", header=F, sep=",", stringsAsFactors =F)
yoshihara.candidates <- read.delim("Data/Candidate.Genes/Yoshihara2012.sig.mapped.csv", header=F, sep=",", stringsAsFactors =F)
candidates <- unique(unlist(c(tcga.candidates, yoshihara.candidates)))

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

png("Figures/90th.Quantile.Gene.Expression.Histogram")
hist(TCGA.quantile, breaks=50, main="Distribution of Expression 90th Quantile Thresholds",
     xlab="90th Quantile Threshold",
     ylab="Number of Genes") 
abline(v=5, col="red", lty="dashed")
dev.off()


TCGA.quantile.mat <- TCGA[TCGA.quantile > 5,]
TCGA.quantile.candidates <- TCGA[TCGA.quantile > 5 | rownames(TCGA) %in% candidates,]
TCGA.quantile.cor <- getCorrelationMatrix(t(TCGA.quantile.mat))
TCGA.quantile.candidates.cor <- getCorrelationMatrix(t(TCGA.quantile.candidates))


for(i in 1:length(args))
{
  TCGA.quantile.cor.bin <- corrMat2BinMat(TCGA.quantile.cor, as.numeric(args[i]))
  fname <- paste("TCGA.Quantile.cor.bin.", args[i], ".txt", sep="")
  write.table(TCGA.quantile.cor.bin, file.path("Data", fname), sep="\t")
  
  TCGA.quantile.candidates.cor.bin <- corrMat2BinMat(TCGA.quantile.candidates.cor, as.numeric(args[i]))
  fname <- paste("TCGA.Quantile.candidates.cor.bin.", args[i], ".txt", sep="")
  write.table(TCGA.quantile.candidates.cor.bin, file.path("Data", fname), sep="\t")
  
  
}

#Iterate over the redundancy values in the parameter sweep
for(i in redundancy)
{
  #Plot the number of eligible genes by correlation threshold
  #Highlight teh correlation thresholds used in the parameter sweep
  rPs <- seq(0.5, 0.99, 0.01)
  nEligible <- getEligibleGenes(TCGA.quantile.cor, rPs, i)
  colors <- rep("grey", length(rPs))
  colors[rPs == 0.6 | rPs == 0.65 | rPs == 0.7] <- "red"
  fname <- paste("Figures/NumberEligibleGenesByRP.Redundancy-", i, ".png",sep="")
  png(fname, width=700, height=900)
  barplot(nEligible, names.arg=rPs, col=colors, main="Number of Eligible Genes by |rP| threshold",
          xlab="|rP|", ylab="Number of Eligible Genes", cex.lab=1.5, cex.main=2)
  dev.off()
}




