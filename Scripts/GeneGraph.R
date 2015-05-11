#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LIBRARIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(igraph)
library(curatedOvarianData)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Get the directly measured genes when DM size = 400, Redundancy = 1, and Correlation Threshold = 0.7
measured <- read.table("Data/Quantile.GGS.Parameter.Sweep.Results/genesets/TCGA.1.400.0.7.measured.txt", sep=",")[1,]

#Get the predictable genes when DM size = 400, Redundancy = 1, and Correlation Threshold = 0.7
predicted <- read.table("Data/Quantile.GGS.Parameter.Sweep.Results/genesets/TCGA.1.400.0.7.covered.txt", sep=",")[1,]

#get the good sample IDs
tcga.samples <- read.csv("Data/SampleFilter/TCGA_eset_samplesRemoved_noDoppel.csv")

#Load the TCGA data
data(TCGA_eset)
TCGA <- exprs(TCGA_eset)
TCGA <- TCGA[,as.character(tcga.samples$x)]
rm(TCGA_eset)
gc()


#Calculate the 90th quantile for each gene in the TCGA data
quantiles <- unlist(apply(TCGA, 1, function(x){
  return(quantile(x, probs=c(0.9)))
}))

#Get gene names for those genes with a 90th quantile threshold >= 5
genes.quant <- rownames(TCGA)[quantiles >= 5]

#Crate a correlation matrix
TCGA.cor <- cor(t(TCGA[genes.quant,]))

#Make the correlation matrix binary
TCGA.cor.bin <- corrMat2BinMat(TCGA.cor, 0.7)

#Zero out the upper triangle of the binary matrix since
#we do not want directed edges.
TCGA.cor.bin[upper.tri(TCGA.cor.bin, diag=T)] <- 0

#Create the graph
g <- graph.adjacency(TCGA.cor.bin)

#Annotate each nodem (gene) to indicate whether it is measured, predicted, or NA
V(g)$Set <- "NA"
V(g)$Set[V(g)$name %in% as.character(unlist(c(measured[1,])))] <- "Measured"
V(g)$Set[V(g)$name %in% as.character(unlist(c(predicted[1,])))] <- "Predicted"

#Indicate the color of each node (gene)
V(g)$color <- ifelse(V(g)$Set == "Measured", "Red", "Blue")

#Save the graph so it can be plotted in Cytoskape
write.graph(g, "Data/TCGA.0.7.GGS.sets.gml", format="gml")




