library("RTCGAToolbox")

#Get the RNAseq and Microarray Data
BrCa.TCGA <- getFirehoseData("BRCA", runDate="20150402", RNAseq2_Gene_Norm = TRUE)


#Get the 95th quantile of expression values for each gene
quantileThresh <- apply(BrCa.TCGA@RNASeq2GeneNorm, 1, function(x){
  thresh <- quantile(x, probs=c(0.95), na.rm=TRUE)
})


trulyExpressedGenes <- quantileThresh > 100

#Create the gene-by-gene correlation matrix for the RNAseq data
BrCa.RNAseq.cor <- cor(t(BrCa.TCGA@RNASeq2GeneNorm[trulyExpressedGenes,]))

#Save the correlation matrix
write.table(BrCa.RNAseq.cor, "BrCa.RNAseq.filtered.cor.txt", sep="\t", quote=FALSE)

filterCorMatrix <- function(inputMatrix, threshold)
{
  m <- nrow(inputMatrix)
  result <- matrix(rep(0,m*m), nrow=m)
  result[inputMatrix >= threshold] <- 1
  colnames(result) <- rownames(result) <- rownames(inputMatrix)
  return(result)
  
}


#Filter the matrix to only include gene-gene correlations > 0.60, 0.65, and 0.70
cor.0.60 <- filterCorMatrix(BrCa.RNAseq.cor, 0.6)
write.table(cor.0.60, "BrCa.RNAseq.filtered.cor.0.60.txt", quote=FALSE, sep="\t")
rm(cor.0.60)
gc()

cor.0.65 <- filterCorMatrix(BrCa.RNAseq.cor, 0.65)
write.table(cor.0.65, "BrCa.RNAseq.filtered.cor.0.65.txt", quote=FALSE, sep="\t")
rm(cor.0.65)
gc()

cor.0.70 <- filterCorMatrix(BrCa.RNAseq.cor, 0.7)
write.table(cor.0.70, "BrCa.RNAseq.filtered.cor.0.70.txt", quote=FALSE, sep="\t")
rm(cor.0.70)
gc()








