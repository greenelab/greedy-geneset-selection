#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LIBRARIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(curatedOvarianData) # >= 1.3.3
library("doppelgangR")
library(reshape2)
library(outliers)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#CONSTANTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vars <- c("sample_type", "histological_type", "grade", "primarysite", "arrayedsite", "summarystage", "tumorstage", "substage", "pltx", "tax", "neo", "recurrence_status", "vital_status", "os_binary", "relapse_binary", "site_of_tumor_first_recurrence", "primary_therapy_outcome_success", "debulking")
minimumSamples <- 120

#Esets that should never be included: Dressman, Bentink, and Pils datasets
excludeEsets <- c("PMID17290060_eset", "E.MTAB.386_eset", "GSE49997_eset")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#FUNCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

getAllDataSets <- function(packageName)
{
  tmp <- data(package=packageName)
  datasets <- tmp[3][[1]][,3]
  return(datasets)
}

exclusionTable <- function(esetVec)
{
  outTable <- data.frame()
  sampleList <- list()
  for(i in 1:length(esetVec))
  {
    exprs <- paste("data(",esetVec[i],")",sep="")
    eval(parse(text=exprs))
    
    exprs <- paste("pheno <- pData(",esetVec[i],")",sep="")
    eval(parse(text=exprs))
    
    goodSamples <- rownames(pheno)
    nTumor <- nSerous <- nEndo <- nHighSerous <- nLowSerous <- nHighEndo <- 0
    
    
    nTumor <- sum(pheno$sample_type == "tumor", na.rm=T)
    goodSamples <- rownames(pheno)[pheno$sample_type == "tumor" %in% TRUE]
    #     nHealthy <- sum(pheno$sample_type == "healthy", na.rm=T)
    #     nBenign <- sum(pheno$sample_type == "benign", na.rm=T)
    #     nMetastatic <- sum(pheno$sample_type == "metastatic", na.rm=T)
    #     nCellline <- sum(pheno$sample_type == "cellline", na.rm=T)
    #     nBorderline <- sum(pheno$sample_type == "borderline", na.rm=T)
    
    if(sum(is.na(pheno$histological_type), na.rm=T) < nrow(pheno))
    {
      #histological type must be assigned
      nSerous <- sum(pheno$sample_type == "tumor" & pheno$histological_type == "ser", na.rm=T)
      nEndo <- sum(pheno$sample_type == "tumor" & pheno$histological_type == "endo", na.rm=T)
      
      goodSamples <- rownames(pheno)[((pheno$sample_type == "tumor" & pheno$histological_type == "ser") | (pheno$sample_type == "tumor" & pheno$histological_type == "endo"))  %in% TRUE]
      nMissingGrade <- sum(pheno$sample_type == "tumor" & pheno$histological_type == "ser" & is.na(pheno$grade), na.rm=T)
      if(sum(is.na(pheno$grade), na.rm=T) < nrow(pheno))
      {
        #grade must be assigned
        nHighSerous <- sum(pheno$sample_type == "tumor" & pheno$histological_type == "ser" & pheno$grade > 1, na.rm=T)
        nLowSerous <- sum(pheno$sample_type == "tumor" & pheno$histological_type == "ser" & pheno$grade == 1, na.rm=T)
        nHighEndo <- sum(pheno$sample_type == "tumor" & pheno$histological_type == "endo" & pheno$grade > 2, na.rm=T)
        goodSamples <- rownames(pheno)[((pheno$sample_type == "tumor" & pheno$histological_type == "ser" & pheno$grade > 1) | (pheno$sample_type == "tumor" & pheno$histological_type == "endo" & pheno$grade > 2)) %in% TRUE]
      }
      
    }
    totSamples <- nrow(pheno)
    #goodSamples <- rownames(pheno)[((pheno$sample_type == "tumor" & pheno$histological_type == "ser" & pheno$grade > 1) | (pheno$sample_type == "tumor" & pheno$histological_type == "endo" & pheno$grade > 2)) %in% TRUE]
    sampleList[[i]] <- NA
    if(length(goodSamples) > 1)
    {
      sampleList[[i]] <- goodSamples
    }
    
    thisInclusion <- c(totSamples, nTumor, nSerous, nMissingGrade, nHighSerous, nLowSerous, nEndo, nHighEndo, length(goodSamples))
    
    if(ncol(outTable) < 1)
    {
      outTable <- data.frame(X=thisInclusion)
      #colnames(outTable) <- esetVec[i]
    }
    
    else
    {
      outTable <- cbind(outTable, thisInclusion)
    }
    
    exprs <- paste("rm(",esetVec[i],", envir=.GlobalEnv)",sep="")
    eval(parse(text=exprs))
  }
  colnames(outTable) <- esetVec
  rownames(outTable) <- c("TotalSamples", "Tumor", "Serous", "Missing Grade","HighGradeSerous", "LowGradeSerous", "Endo", "HighGradeEndo", "AnalyticSet")
  names(sampleList) <- esetVec
  output <- list(outTable, sampleList)
  return(output)
}

simpleExclusion <- function(eset)
{
  
  pheno <- pData(eset)
  
  goodSamples <- rownames(pheno)
  nTumor <- nSerous <- nEndo <- nHighSerous <- nLowSerous <- nHighEndo <- 0
  
  
  nTumor <- sum(pheno$sample_type == "tumor", na.rm=T)
  goodSamples <- rownames(pheno)[pheno$sample_type == "tumor" %in% TRUE]
  #     nHealthy <- sum(pheno$sample_type == "healthy", na.rm=T)
  #     nBenign <- sum(pheno$sample_type == "benign", na.rm=T)
  #     nMetastatic <- sum(pheno$sample_type == "metastatic", na.rm=T)
  #     nCellline <- sum(pheno$sample_type == "cellline", na.rm=T)
  #     nBorderline <- sum(pheno$sample_type == "borderline", na.rm=T)
  
  if(sum(is.na(pheno$histological_type), na.rm=T) < nrow(pheno))
  {
    #histological type must be assigned
    nSerous <- sum(pheno$sample_type == "tumor" & pheno$histological_type == "ser", na.rm=T)
    nEndo <- sum(pheno$sample_type == "tumor" & pheno$histological_type == "endo", na.rm=T)
    
    goodSamples <- rownames(pheno)[((pheno$sample_type == "tumor" & pheno$histological_type == "ser") | (pheno$sample_type == "tumor" & pheno$histological_type == "endo"))  %in% TRUE]
    nMissingGrade <- sum(pheno$sample_type == "tumor" & pheno$histological_type == "ser" & is.na(pheno$grade), na.rm=T)
    if(sum(is.na(pheno$grade), na.rm=T) < nrow(pheno))
    {
      #grade must be assigned
      nHighSerous <- sum(pheno$sample_type == "tumor" & pheno$histological_type == "ser" & pheno$grade > 1, na.rm=T)
      nLowSerous <- sum(pheno$sample_type == "tumor" & pheno$histological_type == "ser" & pheno$grade == 1, na.rm=T)
      nHighEndo <- sum(pheno$sample_type == "tumor" & pheno$histological_type == "endo" & pheno$grade > 2, na.rm=T)
      goodSamples <- rownames(pheno)[((pheno$sample_type == "tumor" & pheno$histological_type == "ser" & pheno$grade > 1) | (pheno$sample_type == "tumor" & pheno$histological_type == "endo" & pheno$grade > 2)) %in% TRUE]
    }
    
  }
  totSamples <- nrow(pheno)
  #goodSamples <- rownames(pheno)[((pheno$sample_type == "tumor" & pheno$histological_type == "ser" & pheno$grade > 1) | (pheno$sample_type == "tumor" & pheno$histological_type == "endo" & pheno$grade > 2)) %in% TRUE]
  sampleList <- NA
  if(length(goodSamples) > 1)
  {
    sampleList <- goodSamples
  }
  
  thisInclusion <- c(totSamples, nTumor, nSerous, nMissingGrade, nHighSerous, nLowSerous, nEndo, nHighEndo, length(goodSamples))
  outTable <- data.frame(X=thisInclusion)
  
  result <- list()
  result[[1]] <- outTable
  result[[2]] <- sampleList
  return(result)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ANALYSIS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#All the datasets within the curatedOvarianData package
esets <- getAllDataSets("curatedOvarianData")

#load the Konecny data
load("Data/Mayo/Mayo.Rda")

#Use the inclusion/exclusion decision tree to filter samples in all curatedOvarainData datasets
inclusionTable <- exclusionTable(esets)

#Use the inclusion/exclusion decision tree on the Konecny data
inclusionTable.mayo <- simpleExclusion(mayo.eset)

#Combine the inclusion results from curatedOvarianData and Konecny
inclusionTable[[1]] <- cbind(inclusionTable[[1]], inclusionTable.mayo[[1]])

#The first list element in the inclusionTable is the data.frame which details the creation of the analytic set and how many samples were excluded and why
colnames(inclusionTable[[1]])[ncol(inclusionTable[[1]])] <- "Mayo.eset"

#Save the data.frame to the harddrive
write.csv(inclusionTable[[1]], "Data/Inclusions.csv")

#The second list element in inclusionTable is a list of dataset specific 'good' sample IDs; i.e. samples which should be included
#JR-Note: remove the following line
#inclusionTable[[2]][[length(inclusionTable[[2]]) + 1]] <- inclusionTable.mayo[[2]]
goodSamples <- inclusionTable[[2]]

#Remove the extra TCGA data (rnaseq and mirna)
goodSamples <- goodSamples[-1 * grep("rna|RNA", names(goodSamples))]

#Only consider the esets with the minimum number of samples
esetList.chosen <- list()
goodSamples.chosen <- list()
for(i in 1:(length(goodSamples)))
{
  if(length(goodSamples[[i]]) >= minimumSamples & !(names(goodSamples)[i] %in% excludeEsets))
  {
    #get the samples
    goodSamples.chosen[[names(goodSamples)[i]]] <- goodSamples[[i]]
    
    #load the eset
    exprs <- paste("data(",names(goodSamples)[i],")",sep="")
    eval(parse(text=exprs))
    
    rm(exprs)
    
    #Limit it to only the good samples
    exprsString <- paste(names(goodSamples)[i], " <- ",names(goodSamples)[i],"[,goodSamples[[i]]]",sep="" )
    eval(parse(text=exprsString))
    
    #     exprsString <- paste("exprs(",names(goodSamples)[i],") <- exprs(", names(goodSamples)[i], ")[,goodSamples[[i]]]",sep="" )
    #     eval(parse(text=exprsString))
    #     exprsString <- paste("pData(",names(goodSamples)[i],") <- pData(", names(goodSamples)[i], ")[goodSamples[[i]],]",sep="" )
    #     eval(parse(text=exprsString))
    
    #add the eset to the eset List
    exprsString <- paste("esetList.chosen[[",length(goodSamples.chosen),"]] <- ",names(goodSamples)[i], sep="")
    eval(parse(text=exprsString))
    
    
    #delete the eset
    exprsString <- paste("rm(",names(goodSamples)[i],")",sep="")
    eval(parse(text=exprsString))
  }
}
names(esetList.chosen) <- names(goodSamples.chosen)

#Add the konecny data to the esetList.chosen and goodSamples.chosen
esetList.chosen[[length(esetList.chosen)+1]] <- mayo.eset[,inclusionTable.mayo[[2]]]
goodSamples.chosen[[length(esetList.chosen)]] <- inclusionTable.mayo[[2]]
names(esetList.chosen)[length(esetList.chosen)] <- names(goodSamples.chosen)[length(esetList.chosen)] <- "mayo.eset"

#Pre-process teh esets to improve matching
testesets <- lapply(esetList.chosen, function(X){  
  X$alt_sample_name <- paste(X$sample_type, gsub("[^0-9]", "", X$alt_sample_name), sep="_")
  pData(X) <- pData(X)[, !grepl("uncurated_author_metadata", colnames(pData(X)))]
  return(X) })

#Use doppelgangR to find similar sample pairs accross datasets except for mayo; this code assumes that the mayo data is the last eset in the list
result1 <- doppelgangR(testesets[-1 * length(testesets)], corFinder.args=list(use.ComBat=TRUE), cache.dir=NULL)

#Process the doppelgangR results into data.frames and write to the harddrive
result1.df <- summary(result1)
result1.full <- result1@fullresults
result1.df.full <- c()
for(i in 1:length(result1.full))
{
  tmp <- merge(result1.full[[i]][["expr.doppels"]][["outlierFinder.res"]], result1.full[[i]][["pheno.doppels"]][["outlierFinder.res"]], by=c("sample1", "sample2"))
  colnames(tmp) <- c("sample1", "sample2", "expr.similarity", "expr.doppel", "pheno.similarity", "pheno.doppel")
  result1.df.full <- rbind(result1.df.full, tmp)
}
write.table(result1.df.full, file="Data/doppelgangR/pairwiseSampleComparisons.tsv", sep="\t", quote=FALSE, row.names=FALSE)


#Use doppelgangR to find similar sample pairs in TCGA and konecny data sets using only the expression data
data(TCGA_eset)
result2 <- doppelgangR(list(TCGA_eset=TCGA_eset, mayo.eset=mayo.eset), corFinder.args=list(use.ComBat=TRUE), phenoFinder.args=NULL, cache.dir=NULL)

#Process the doppelgangR results into data.frames and write to the harddrive
result2 <- doppelgangR(list(TCGA_eset=TCGA_eset[,inclusionTable[[2]]["TCGA_eset"]], mayo.eset=mayo.eset[,inclusionTable.mayo[[2]]]), corFinder.args=list(use.ComBat=TRUE), phenoFinder.args=NULL, cache.dir=NULL)
result2.df <- summary(result2)
result2.full <- result2@fullresults
result2.df.full <- c()
for(i in 1:length(result2.full))
{
  tmp <- merge(result2.full[[i]][["expr.doppels"]][["outlierFinder.res"]], result2.full[[i]][["pheno.doppels"]][["outlierFinder.res"]], by=c("sample1", "sample2"))
  colnames(tmp) <- c("sample1", "sample2", "expr.similarity", "expr.doppel", "pheno.similarity", "pheno.doppel")
  result2.df.full <- rbind(result2.df.full, tmp)
}
write.table(result2.df.full, file="Data/doppelgangR/pairwiseSampleComparisons_MayoTCGA.tsv", sep="\t", quote=FALSE, row.names=FALSE)

#Comebine the doppelgangR results and write to harddrive
tmp <- rbind(result1.df.full, result2.df.full)
write.table(tmp[!duplicated(tmp[,c(1,2)]),], file="Data/doppelgangR/pairwiseSampleComparisons_Combined.tsv", sep="\t", quote=FALSE, row.names=FALSE)



#Get list of expression doppelgangR samples; we assume that if they are doppelgangR sample pairs and they are more than 0.95 simlar then they are duplicates
doppelSamples <- c(result1.df$sample1[result1.df$expr.doppel & result1.df$expr.similarity > 0.95], result1.df$sample2[result1.df$expr.doppel & result1.df$expr.similarity > 0.95],
                   result2.df$sample1[result2.df$expr.doppel & result2.df$expr.similarity > 0.95], result2.df$sample2[result2.df$expr.doppel & result2.df$expr.similarity > 0.95])
#get rid of the X1: prefixes
doppelSamples <- unlist(strsplit(doppelSamples,":"))[seq(2,length(doppelSamples)*2,2)]
doppelSamples <- unique(doppelSamples)

#write out the 'good' samples
for(i in 1:length(goodSamples.chosen))
{
  sampleList <- setdiff(goodSamples.chosen[[i]], doppelSamples)
  outFName <- paste("Data/SampleFilter/", names(goodSamples.chosen)[i], "_samplesRemoved_noDoppel.csv", sep="")
  write.csv(sampleList, outFName)
  cat(outFName)
  cat("\n")
}


goodSamples.doppelRemoved <- setdiff(as.vector(unlist(goodSamples.chosen)), doppelSamples)
numberRemoved <- length(as.vector(unlist(goodSamples.chosen))) - length(goodSamples.doppelRemoved)
print(paste("Removed", numberRemoved,"samples based on expression similarity."))
minimumExpressionSimilarity <- min(result1.df$expr.similarity[result1.df$expr.doppel & result1.df$expr.similarity > 0.95])
print(paste("Minimum expression similarity in duplicates:", minimumExpressionSimilarity))