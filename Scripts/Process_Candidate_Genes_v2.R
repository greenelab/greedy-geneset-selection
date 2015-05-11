
library(curatedOvarianData)
library(HGNChelper)

data(TCGA_eset)
TCGA <- exprs(TCGA_eset)

#Yoshihara candidate genes
yoshihara.genes <- read.delim("Data/Candidate.Genes/Yoshihara2012.sig.csv", sep=",")
yoshihara.genes.corrections <- checkGeneSymbols(yoshihara.genes$Gene.Symbol)
yoshihara.genes.corrected <- yoshihara.genes.corrections$Suggested.Symbol


sum(rownames(TCGA) %in% yoshihara.genes.corrected)
yoshihara.genes.final <- rownames(TCGA)[rownames(TCGA) %in% yoshihara.genes.corrected]

#TCGA candidate genes
tcga.genes <- read.delim("Data/Candidate.Genes/TCGA2011.sig.csv", sep=",")
tcga.genes.corrections <- checkGeneSymbols(tcga.genes$Name)
tcga.genes.corrected <- tcga.genes.corrections$Suggested.Symbol


sum(rownames(TCGA) %in% tcga.genes.corrected)
tcga.genes.final <- rownames(TCGA)[rownames(TCGA) %in% tcga.genes.corrected]

write.table(yoshihara.genes.final, "Data/Candidate.Genes/Yoshihara2012.sig.mapped.csv", sep=",", quote=F, row.names=F, col.names=F, eol=",")
write.table(tcga.genes.final, "Data/Candidate.Genes/TCGA2011.sig.mapped.csv", sep=",", quote=F, row.names=F, col.names=F, eol=",")