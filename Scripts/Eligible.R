#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LIBRARIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


library(ff)
library(doMC)
library(curatedOvarianData)
library(wq)
library(ggplot2)
library(scales)

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


TCGA <- TCGA[TCGA.quantile > 5,]
TCGA.cor <- getCorrelationMatrix(t(TCGA))


plots <- list()
thresh <- seq(0.5, 1, 0.01)
colors <- rep("lightgrey", length(thresh))
colors[thresh %in% c(0.6, 0.65, 0.7)] <- "red"
for(i in 1:3)
{
  
  eg <- getEligibleGenes(TCGA.cor, thresh, i)
  plots[[i]] <- data.frame(thresh, eg)
  print(paste("Finished ", i, sep=" "))
}

p1 <- ggplot(plots[[1]], aes(x=thresh, y=eg)) + geom_bar(stat="identity", col="black", fill=colors) + theme_bw() +
  xlab(bquote(abs(r[p]))) +
  ylab("Number of Eligible Genes") +
  scale_y_continuous(limits=c(0,6100), oob=rescale_none, breaks=seq(0,6000,1000)) +
  theme(text = element_text(size=20)) +
  labs(title="Min. correlated genes=1")

p2 <- ggplot(plots[[2]], aes(x=thresh, y=eg)) + geom_bar(stat="identity", col="black", fill=colors) + theme_bw() +
  xlab(bquote(abs(r[p]))) +
  ylab("Number of Eligible Genes") +
  scale_y_continuous(limits=c(0,6100), oob=rescale_none, breaks=seq(0,6000,1000)) +
  theme(text = element_text(size=20)) +
  labs(title="Min. correlated genes=2")

p3 <- ggplot(plots[[3]], aes(x=thresh, y=eg)) + geom_bar(stat="identity", col="black", fill=colors) + theme_bw() +
  xlab(bquote(abs(r[p]))) +
  ylab("Number of Eligible Genes") +
  scale_y_continuous(limits=c(0,6100), oob=rescale_none, breaks=seq(0,6000,1000)) +
  theme(text = element_text(size=20)) +
  labs(title="Min. correlated genes=3")

# layOut(list(p1, 1, 1),
#        list(p2, 1, 2),
#        list(p3,1,3),
#        list(plot.quantile,2,1:3))

plot.quantile2 <- ggplot(data.frame(TCGA.quantile), aes(x=TCGA.quantile)) + geom_density(colour="black", fill="blue", size=2) +  theme_bw() +
  scale_x_continuous(limits=c(0,round(max(TCGA.quantile)) + 1), breaks=seq(0,round(max(TCGA.quantile)) + 1,2)) +
  geom_vline(aes(xintercept=5), color="red", linetype="dashed", size=2) +
  xlab("90th Quantile Threshold") +
  ylab("Frequency") +
  theme(text = element_text(size=20))

png("Figures/Eligible.Genes.Plots.png", width=1400,height=790)
layOut(list(p1, 2, 1),
       list(p2, 2, 2),
       list(p3,2,3),
       list(plot.quantile2,1,1:3))
dev.off()







#Genes for enrichment
TCGA.cor.bin <- corrMat2BinMat(TCGA.cor, 0.60)
genes <- colnames(TCGA.cor.bin)[colSums(TCGA.cor.bin) > 1]
genes.split <- unlist(strsplit(genes, "///"))
write.table(genes.split, "genelist.quantile.filtered.0.6.1.txt", quote=F, col.names=F, row.names=F, sep="\n")


#Read in the GO-slimo pathways 
pathways <- read.csv2("Data/OBO.genelist.csv", header=F, sep=";")

#Background genes
enrichment.df <- data.frame(GO_ID=pathways[,1], Pathway=pathways[,2])

enrichment.df$backgroundFreq <- NA
enrichment.df$sampleFreq <- NA
enrichment.df$Expected <- NA
enrichment.df$pvalue <- NA

genes.background <- unlist(strsplit(colnames(TCGA.cor.bin), "///"))
genes.eligible <- genes.split
for(i in 1:nrow(pathways))
{
  genes.pathway <- gsub(" ", "", unlist(strsplit(as.character(pathways[i,3]), ",")))
  eligible.overlap <- length(intersect(genes.eligible, genes.pathway))
  background.overlap <- length(intersect(genes.background, genes.pathway))
  nAllGenes <- length(genes.background)
  nEligibleGenes <- length(genes.eligible)
  
  enrichment.df$backgroundFreq[i] <- background.overlap
  enrichment.df$sampleFreq[i] <- eligible.overlap
  
  expected <- (length(genes.eligible) / length(genes.background)) * background.overlap
  enrichment.df$Expected[i] <- round(expected)

  p <- phyper(eligible.overlap, length(genes.eligible), length(genes.background) - length(genes.eligible), background.overlap, lower.tail=FALSE )
  enrichment.df$pvalue[i] <- as.character(p)
  
}

enrichment.df$pAdjusted <- p.adjust(enrichment.df$pvalue, method="bonferroni")
enrichment.df <- enrichment.df[order(enrichment.df$pAdjusted, decreasing=F),]

write.table(enrichment.df, "Data/TCGA.1.0.6.eligible.GoSlim.enrichment.table.csv", quote=F, sep=",", row.names=F)
