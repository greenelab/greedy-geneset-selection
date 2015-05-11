library(HGNChelper)
library(curatedOvarianData)




candidate.genes.TCGA <- read.csv("~/Documents/tmp/evolutionary-gene-set-selection/Data/Candidate.Genes/TCGA2011.sig.csv", as.is=T)

data(TCGA_eset)
TCGA <- exprs(TCGA_eset)
rm(TCGA_eset)

hgnc.corrections <- checkGeneSymbols(candidate.genes.TCGA$Name)
hgnc.corrections[!hgnc.corrections$Approved,]

candidate.genes.TCGA$Name.Corrected <- hgnc.corrections$Suggested.Symbol

TCGA.candidate.genes <- intersect(rownames(TCGA), candidate.genes.TCGA$Name.Corrected)
tmp <- c()
tmp <- rbind(tmp, TCGA.candidate.genes)
write.table(tmp, "~/Documents/tmp/evolutionary-gene-set-selection/Data/Candidate.Genes/TCGA2011.sig.processed.csv", quote=F, row.names=F, col.names=F, sep=",")
rm(tmp, TCGA.candidate.genes, hgnc.corrections)
gc()







candidate.genes.Yoshihara <- read.csv("~/Documents/tmp/evolutionary-gene-set-selection/Data/Candidate.Genes/Yoshihara2012.sig.csv", as.is=T)
hgnc.corrections <- checkGeneSymbols(candidate.genes.Yoshihara$Gene.Symbol)
hgnc.corrections[!hgnc.corrections$Approved,]

candidate.genes.Yoshihara$Name.Corrected <- hgnc.corrections$Suggested.Symbol
Yoshihara.candidate.genes <- intersect(rownames(TCGA), candidate.genes.Yoshihara$Name.Corrected)
tmp <- c()
tmp <- rbind(tmp, Yoshihara.candidate.genes)
write.table(tmp, "~/Documents/tmp/evolutionary-gene-set-selection/Data/Candidate.Genes/Yoshihara2012.sig.processed.csv", quote=F, row.names=F, col.names=F, sep=",")

