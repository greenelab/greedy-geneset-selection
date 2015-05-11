

TCGA.sig <- read.csv("Data/Candidate.Genes/TCGA2011.sig.csv", as.is=T)
symbol.map <- checkGeneSymbols(TCGA.sig$Name) #using HGNChelper
symbol.map[!symbol.map$Approved,]
TCGA.sig$Names.suggested <- symbol.map$Suggested.Symbol


Yoshi2012.sig <- read.csv("Data/Candidate.Genes/Yoshi2012.sig.csv", as.is=T)
hgnc.corrections <- checkGeneSymbols(Yoshi2012.sig$Gene.Symbol)
hgnc.corrections[!hgnc.corrections$Approved,]
Yoshi2012.sig$Corrected.Gene.Symbol <- hgnc.corrections$Suggested.Symbol
