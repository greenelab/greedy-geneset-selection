args <- commandArgs(trailingOnly=TRUE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LIBRARIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ggplot2)
library(grid)
library(curatedOvarianData)

tcga.candidates <- read.delim("Data/Candidate.Genes/TCGA2011.sig.mapped.csv", header=F, sep=",", stringsAsFactors =F)
yoshihara.candidates <- read.delim("Data/Candidate.Genes/Yoshihara2012.sig.mapped.csv", header=F, sep=",", stringsAsFactors =F)
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

process_results <- function(fName, binLevel, plotName)
{
	GGS.results <- read.table(fname, sep="\t", header=T)
	GGS.results$Redundancy <- factor(GGS.results$Fold)

	fname <- paste("Figures/",plotName, binLevel, ".png", sep="")
# 	png(fname)
	ggplot(GGS.results, aes(x=nMeasured, y=nCovered, group=Redundancy,color=Redundancy)) + geom_line() + geom_point() +
  		xlab("# Directly Measured Genes") + ylab("# Predictable Genes")
# 	dev.off()

  ggsave(fname)
	GGS.results$Threshold <- binLevel
	return(GGS.results)
}

threshold_labeler <- function(variable, value)
{
  if(value == "0.6")
  {
    return("|rP| = 0.60")
  }
  
  if(value == "0.66")
  {
    return("|rP| = 0.66")
  }
  
  if(value == "0.7")
  {
    return("|rP| = 0.70")
  }

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
#ANALYSIS
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



TCGA.quantile.mat <- TCGA[TCGA.quantile > 5,]
TCGA.quantile.yoshihara.candidates <- TCGA[TCGA.quantile > 5 | rownames(TCGA) %in% yoshihara.candidates,]
TCGA.quantile.tcga.candidates <- TCGA[TCGA.quantile > 5 | rownames(TCGA) %in% tcga.candidates,]
TCGA.quantile.cor <- getCorrelationMatrix(t(TCGA.quantile.mat))
TCGA.quantile.yoshihara.candidates.cor <- getCorrelationMatrix(t(TCGA.quantile.yoshihara.candidates))
TCGA.quantile.tcga.candidates.cor <- getCorrelationMatrix(t(TCGA.quantile.tcga.candidates))

results <- list()
results.yoshihara <- list()
results.tcga <- list()

for(i in 1:length(args))
{

	#Process and plot the no-candidate gene set results
	fname <- paste("Data/Quantile.GGS.Parameter.Sweep.Results/process_results.", args[i], ".results.txt", sep="")
	GGS.results <- process_results(fname, args[i], "CoverageByNumGenesMeasured.TCGA.Quantile.")
	results[[i]] <- GGS.results

	#Process and plot the Yoshihara candidate gene set results
	fname <- paste("Data/Quantile.GGS.Parameter.Sweep.Results/Yoshihara.Candidate.Genes/process_results.", args[i], ".results.txt", sep="")
	GGS.results <- process_results(fname, args[i], "CoverageByNumGenesMeasured.TCGA.Quantile.Yoshihara.Candidate.Genes.")
	results.yoshihara[[i]] <- GGS.results


	#Process and plot the TCGA candidate gene set results
	fname <- paste("Data/Quantile.GGS.Parameter.Sweep.Results/TCGA.Candidate.Genes/process_results.", args[i], ".results.txt", sep="")
	GGS.results <- process_results(fname, args[i], "CoverageByNumGenesMeasured.TCGA.Quantile.TCGA.Candidate.Genes.")
	results.tcga[[i]] <- GGS.results

}

results.combined <- do.call("rbind", results)
results.yoshihara.combined <- do.call("rbind", results.yoshihara)
results.tcga.combined <- do.call("rbind", results.tcga)

results.combined$eligible[results.combined$Threshold == "0.6" & results.combined$Fold == 1] <- getEligibleGenes(TCGA.quantile.cor, 0.6, 1)
results.combined$eligible[results.combined$Threshold == "0.6" & results.combined$Fold == 2] <- getEligibleGenes(TCGA.quantile.cor, 0.6, 2)
results.combined$eligible[results.combined$Threshold == "0.6" & results.combined$Fold == 3] <- getEligibleGenes(TCGA.quantile.cor, 0.6, 3)

results.combined$eligible[results.combined$Threshold == "0.65" & results.combined$Fold == 1] <- getEligibleGenes(TCGA.quantile.cor, 0.65, 1)
results.combined$eligible[results.combined$Threshold == "0.65" & results.combined$Fold == 2] <- getEligibleGenes(TCGA.quantile.cor, 0.65, 2)
results.combined$eligible[results.combined$Threshold == "0.65" & results.combined$Fold == 3] <- getEligibleGenes(TCGA.quantile.cor, 0.65, 3)

results.combined$eligible[results.combined$Threshold == "0.7" & results.combined$Fold == 1] <- getEligibleGenes(TCGA.quantile.cor, 0.7, 1)
results.combined$eligible[results.combined$Threshold == "0.7" & results.combined$Fold == 2] <- getEligibleGenes(TCGA.quantile.cor, 0.7, 2)
results.combined$eligible[results.combined$Threshold == "0.7" & results.combined$Fold == 3] <- getEligibleGenes(TCGA.quantile.cor, 0.7, 3)

results.combined$missed <- pmax(0, results.combined$eligible - (results.combined$nMeasured + results.combined$nCovered))
results.combined$Redundancy <- factor(results.combined$Fold)

fig <- ggplot(results.combined, aes(x=nMeasured, y=nCovered, color=Redundancy)) + theme_bw() + 
  geom_line() + geom_point(aes(shape=Redundancy)) +
  geom_line(aes(x=nMeasured, y=missed, color=Redundancy), linetype="dotted", size=1) + 
  facet_wrap(~ Threshold) + xlab("# Directly Measured Genes") + ylab("# Predictable Genes") +
  theme(text = element_text(size=20), panel.margin=unit(1.5, "lines")) 

ggsave("Figures/CoverageByNumGenesMeasured.Quantile.png", fig, width=11, height=6)


results.yoshihara.combined$eligible[results.yoshihara.combined$Threshold == "0.6" & results.yoshihara.combined$Fold == 1] <- getEligibleGenes(TCGA.quantile.yoshihara.candidates.cor, 0.6, 1)
results.yoshihara.combined$eligible[results.yoshihara.combined$Threshold == "0.6" & results.yoshihara.combined$Fold == 2] <- getEligibleGenes(TCGA.quantile.yoshihara.candidates.cor, 0.6, 2)
results.yoshihara.combined$eligible[results.yoshihara.combined$Threshold == "0.6" & results.yoshihara.combined$Fold == 3] <- getEligibleGenes(TCGA.quantile.yoshihara.candidates.cor, 0.6, 3)

results.yoshihara.combined$eligible[results.yoshihara.combined$Threshold == "0.65" & results.yoshihara.combined$Fold == 1] <- getEligibleGenes(TCGA.quantile.yoshihara.candidates.cor, 0.65, 1)
results.yoshihara.combined$eligible[results.yoshihara.combined$Threshold == "0.65" & results.yoshihara.combined$Fold == 2] <- getEligibleGenes(TCGA.quantile.yoshihara.candidates.cor, 0.65, 2)
results.yoshihara.combined$eligible[results.yoshihara.combined$Threshold == "0.65" & results.yoshihara.combined$Fold == 3] <- getEligibleGenes(TCGA.quantile.yoshihara.candidates.cor, 0.65, 3)

results.yoshihara.combined$eligible[results.yoshihara.combined$Threshold == "0.7" & results.yoshihara.combined$Fold == 1] <- getEligibleGenes(TCGA.quantile.yoshihara.candidates.cor, 0.7, 1)
results.yoshihara.combined$eligible[results.yoshihara.combined$Threshold == "0.7" & results.yoshihara.combined$Fold == 2] <- getEligibleGenes(TCGA.quantile.yoshihara.candidates.cor, 0.7, 2)
results.yoshihara.combined$eligible[results.yoshihara.combined$Threshold == "0.7" & results.yoshihara.combined$Fold == 3] <- getEligibleGenes(TCGA.quantile.yoshihara.candidates.cor, 0.7, 3)

results.yoshihara.combined$missed <- pmax(0, results.yoshihara.combined$eligible - (results.yoshihara.combined$nMeasured + results.yoshihara.combined$nCovered))
results.yoshihara.combined$Redundancy <- factor(results.yoshihara.combined$Fold)
# fig <- ggplot(results.yoshihara.combined, aes(x=nMeasured, y=nCovered, color=Redundancy)) + theme_bw() +
#   geom_line() + geom_point(aes(shape=Redundancy)) +
#   geom_line(aes(x=nMeasured, y=missed, color=Redundancy), linetype="dotted", size=1) + 
#   facet_wrap(~ Threshold) + xlab("# Directly Measured Genes") + ylab("# Predictable Genes") +
#   theme(text = element_text(size=20))
# ggsave("Figures/CoverageByNumGenesMeasured.Quantile.Yoshihara.Candidate.Genes.png", fig, width=11, height=6)


results.tcga.combined$eligible[results.tcga.combined$Threshold == "0.6" & results.tcga.combined$Fold == 1] <- getEligibleGenes(TCGA.quantile.tcga.candidates.cor, 0.6, 1)
results.tcga.combined$eligible[results.tcga.combined$Threshold == "0.6" & results.tcga.combined$Fold == 2] <- getEligibleGenes(TCGA.quantile.tcga.candidates.cor, 0.6, 2)
results.tcga.combined$eligible[results.tcga.combined$Threshold == "0.6" & results.tcga.combined$Fold == 3] <- getEligibleGenes(TCGA.quantile.tcga.candidates.cor, 0.6, 3)

results.tcga.combined$eligible[results.tcga.combined$Threshold == "0.65" & results.tcga.combined$Fold == 1] <- getEligibleGenes(TCGA.quantile.tcga.candidates.cor, 0.65, 1)
results.tcga.combined$eligible[results.tcga.combined$Threshold == "0.65" & results.tcga.combined$Fold == 2] <- getEligibleGenes(TCGA.quantile.tcga.candidates.cor, 0.65, 2)
results.tcga.combined$eligible[results.tcga.combined$Threshold == "0.65" & results.tcga.combined$Fold == 3] <- getEligibleGenes(TCGA.quantile.tcga.candidates.cor, 0.65, 3)

results.tcga.combined$eligible[results.tcga.combined$Threshold == "0.7" & results.tcga.combined$Fold == 1] <- getEligibleGenes(TCGA.quantile.tcga.candidates.cor, 0.7, 1)
results.tcga.combined$eligible[results.tcga.combined$Threshold == "0.7" & results.tcga.combined$Fold == 2] <- getEligibleGenes(TCGA.quantile.tcga.candidates.cor, 0.7, 2)
results.tcga.combined$eligible[results.tcga.combined$Threshold == "0.7" & results.tcga.combined$Fold == 3] <- getEligibleGenes(TCGA.quantile.tcga.candidates.cor, 0.7, 3)

results.tcga.combined$missed <- pmax(results.tcga.combined$eligible - (results.tcga.combined$nMeasured + results.tcga.combined$nCovered), 0)
results.tcga.combined$Redundancy <- factor(results.tcga.combined$Fold)

# fig <- ggplot(results.tcga.combined, aes(x=nMeasured, y=nCovered, group=Fold,color=Redundancy)) + theme_bw() +
#   geom_line() + geom_point(aes(shape=Redundancy)) +
#   geom_line(aes(x=nMeasured, y=missed, color=Redundancy), linetype="dotted", size=1) + 
#   facet_wrap(~ Threshold) + xlab("# Directly Measured Genes") + ylab("# Predictable Genes") +
#   theme(text = element_text(size=20))
# ggsave("Figures/CoverageByNumGenesMeasured.Quantile.TCGA.Candidate.Genes.png", fig, width=11, height=6)



results.candidates.combined <- rbind(results.tcga.combined, results.yoshihara.combined)
results.candidates.combined$Candidates <- "Yoshihara"
results.candidates.combined$Candidates[1:nrow(results.tcga.combined)] <- "TCGA"
results.candidates.combined$CandidatesII <- factor(results.candidates.combined$Candidates, levels=c("Yoshihara", "TCGA"))

fig <- ggplot(results.candidates.combined, aes(x=nMeasured, y=nCovered, group=Fold,color=Redundancy)) + theme_bw() +
  geom_line() + geom_point(aes(shape=Redundancy)) +
  geom_line(aes(x=nMeasured, y=missed, color=Redundancy), linetype="dotted", size=1) + 
  facet_grid(CandidatesII ~ Threshold) + xlab("# Directly Measured Genes") + ylab("# Predictable Genes") +
  xlim(100,400) + theme(text = element_text(size=20))
ggsave("Figures/CoverageByNumGenesMeasured.Quantile.All.Candidate.Genes.png", fig, width=11, height=11)
