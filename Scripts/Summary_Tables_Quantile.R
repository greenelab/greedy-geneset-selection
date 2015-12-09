#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LIBRARIES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(boot)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#CONSTANTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
allGenesCoverageDir <- "Data/GGS.Parameter.Sweep.Results/"
MADGenesCoverageDir <- "Data/MAD.GGS.Parameter.Sweep.Results/"
QuantileGenesCoverageDir <- "Data/Quantile.GGS.Parameter.Sweep.Results/"

candidateSuffix <- ".Candidate.Genes/"
candidates <- c("Yoshihara", "TCGA")
imputationDir <- "Imputation/Data/"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#FUNCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#http://stackoverflow.com/questions/18341569/r-calculate-the-standard-error-using-bootstrap
meanFunc <- function(x,i)
{
  return(mean(x[i], na.rm=T))
}



subsetImputation <- function(dataOrig, directlyMeasured, correlThresh, candidates, filtering, redundancy)
{
  
  data <- dataOrig[dataOrig$NumberMeasuredGenes == directlyMeasured & 
                     dataOrig$CorrelationThreshold == correlThresh &
                     dataOrig$Candidates == candidates &
                     dataOrig$Filtering == filtering &
                     dataOrig$Redundancy == redundancy,]
  if(nrow(data) < 1)
  {
    res <- list()
    cols <- grep("Spearman", colnames(dataOrig))
    for( i in 1:length(cols))
    {
      name <- colnames(data)[cols[i]]
      res[name] <- NA
    }
    res["TrueNumberPredictedGenes"] <- NA
    return(res)
  }
  else
  {
    res <- list()
    cols <- grep("Spearman", colnames(data))
    for( i in 1:length(cols))
    {
      name <- colnames(data)[cols[i]]
      tmpDat <- data[,cols[i]]
      res[name] <- NA
      if(sum(!is.na(length(tmpDat))) > 0)
      {
        tmpDat.avg <- mean(tmpDat, na.rm=T)
        tmpDat.bt <- boot(tmpDat, meanFunc, 100)
        tmpDat.se <- sd(tmpDat.bt$t)
        res[name] <- paste(round(tmpDat.avg,3)," (", round(tmpDat.se,3), ")",sep="")
        missing.name <- gsub("Spearman", "proportion\\.predicted", name)
        res[missing.name] <- sum(!is.na(tmpDat)) / length(tmpDat)
      }
      
    }
    
    res["TrueNumberPredictedGenes"] <- nrow(data)
    return(res)
  }
  
}



getSummary <- function(covData, impData)
{
  #covData is a list of data frames each summarizing the number of predictable genes as a function of DM size, redundancy, and correlation thresholds
  #        levels of covData are specific to a correlation threshold and the use of filtering and candidate gene sets
  #impData is a data.frame giving the prediction accuracy for each predicted gene
  
  summaryTable <- c()
  for(i in 1:length(covData))
  {
#     corT <- names(covData)[i]
    
    thisDat <- covData[[i]]
    for(j in 1:nrow(thisDat))
    {
      #The FakeNMeasured corresponds to the DM sizes specified in the parameter sweep and is used to create the result file names
      #when using candidate gene sets alone, the size of DM may differ and will be captured by the nMeasured variable
      dirMeas <- thisDat$FakeNMeasured[j]
      trueDirMeas <- thisDat$nMeasured[j]
      #fold <- thisDat$FoldRedundancy[j]
      fold <- thisDat$Fold[j]
      rawScore <- thisDat$Score[j]
      covered <- thisDat$nPredictable[j]
      
      filter <- thisDat$Filtering[j]
      candidate <- thisDat$Candidates[j]
      corT <- thisDat$CorThresh[j]
      
      cat(paste(corT, candidate, filter, fold,trueDirMeas,"       \r",sep="   "))
      predAcc <- subsetImputation(impData, dirMeas, corT, candidate, filter, fold)
      predAcc["CorrelationThreshold"] <- corT
      predAcc["NumberDirectlyMeasured"] <- trueDirMeas
      predAcc["Redundancy"] <- fold
      predAcc["NumberPredicted"] <- covered
      predAcc["Filtering"] <- filter
      predAcc["Candidates"] <- candidate
      #predAcc[["TrueNumberMeasuredGenes"]] <- trueDirMeas


      summaryTable <- rbind(summaryTable, unlist(predAcc))
    }
    
  }
  cat("                                                         \n")
  
  summaryTable <- data.frame(summaryTable)

  summaryTable <- summaryTable[order(summaryTable$Filtering, summaryTable$Candidates, as.numeric(paste(summaryTable$CorrelationThreshold)), as.numeric(paste(summaryTable$Redundancy)), as.numeric(paste(summaryTable$NumberDirectlyMeasured))),]
  #summaryTable <- summaryTable[,c(12,13,8,10,9,11,1,2,3,7,4,5,6)]
  #summaryTableColumnOrder <- c("Filtering", "Candidates", "CorrelationThreshold", "Redundancy", "NumberDirectlyMeasured", "TrueNumberMeasuredGenes", "NumberPredicted", "TrueNumberPredictedGenes", "TCGA.Spearman", "Tothill.Spearman", "Yoshihara.Spearman", "Goode.Spearman", "Bonome.Spearman")
  cols <- colnames(summaryTable)[grep("Spearman", colnames(summaryTable))]
  cols.prop <- colnames(summaryTable)[grep("proportion", colnames(summaryTable))]
  #summaryTableColumnOrder <- c("Filtering", "Candidates", "CorrelationThreshold", "Redundancy", "NumberDirectlyMeasured", "NumberPredicted", "TrueNumberPredictedGenes", cols, cols.prop)
  summaryTableColumnOrder <- c("Filtering", "Candidates", "CorrelationThreshold", "Redundancy", "NumberDirectlyMeasured", "TrueNumberPredictedGenes", cols, cols.prop)
  summaryTable <- summaryTable[,summaryTableColumnOrder]
  
  #print(names(predAcc))
  return(summaryTable)
}

#http://r.789695.n4.nabble.com/Rounding-to-the-nearest-5-td863189.html
mround <- function(x, base)
{
  base*round(x/base)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LOAD DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Let's get the coverage results for Quantile Genes 
files <- list.files(path=QuantileGenesCoverageDir, pattern="process_results.*results.txt", full.names=T, recursive=FALSE)
for(i in 1:length(candidates))
{

  tmp <- list.files(path=paste(QuantileGenesCoverageDir, candidates[i], candidateSuffix,sep=""), pattern="process_results.*results.txt", full.names=T, recursive=FALSE)
  files <- c(files, tmp)

}



#This list is indexed by correlation threshold
coverageResults <- list()
for(i in 1:length(files))
{
  trash <- unlist(strsplit(files[i], "/"))
  fname <- trash[length(trash)]
  params <- unlist(strsplit(fname, "\\."))
  
  if(length(trash) == 4 & trash[length(trash)-2] == unlist(strsplit(allGenesCoverageDir,"/"))[-1])
  {
    filtering = "None"
    thisCand = "None"
  }
  else if(length(trash) == 4 & trash[length(trash)-2] == unlist(strsplit(MADGenesCoverageDir,"/"))[-1])
  {
    filtering = "MAD"
    thisCand = "None"
  }
  else if(length(trash) == 4 & trash[length(trash)-2] == unlist(strsplit(QuantileGenesCoverageDir,"/"))[-1])
  {
    filtering = "Quantile"
    thisCand = "None"
  }
  else if(length(trash) == 5 & trash[length(trash)-3] == unlist(strsplit(allGenesCoverageDir,"/"))[-1])
  {
    filtering = "None"
    for(j in 1:length(candidates))
    {
      if(paste(trash[length(trash)-2],"/",sep="") == paste(candidates[j],candidateSuffix,sep=""))
      {
        thisCand = candidates[j]
      }
    }
  }
  else if(length(trash) == 5 & trash[length(trash)-3] == unlist(strsplit(MADGenesCoverageDir,"/"))[-1])
  {
    filtering = "MAD"
    for(j in 1:length(candidates))
    {
      if(paste(trash[length(trash)-2],"/",sep="") == paste(candidates[j],candidateSuffix,sep=""))
      {
        thisCand = candidates[j]
      }
    }
  }
  else if(length(trash) == 5 & trash[length(trash)-3] == unlist(strsplit(QuantileGenesCoverageDir,"/"))[-1])
  {
    filtering = "Quantile"
    for(j in 1:length(candidates))
    {
      if(paste(trash[length(trash)-2],"/",sep="") == paste(candidates[j],candidateSuffix,sep=""))
      {
        thisCand = candidates[j]
      }
    }
  }
  
  corThresh <- params[3]
  corThresh <- paste("0.",corThresh,sep="")
  
  
  tmp <- read.table(files[i], sep="\t", header=T)
  tmp$Filtering <- filtering
  tmp$Candidates <- thisCand
  tmp$CorThresh <- corThresh
  coverageResults[[i]] <- tmp
  
}





#Let's get the imputation results for All Genes No Candidate Gene Sets
files <- list.files(path=imputationDir, pattern=".*.Quantile.imputation.results.txt$", full.names=T, recursive=FALSE)

imputationResults <- lapply(files, function(x) {
  tmp <- try(read.table(x))
  if(class(tmp) == "try-error")
  {
    return(NA)
  }
  fname <- tail(unlist(strsplit(x, split="/")), n=1)
  trueNGenes <- length(readLines(x)) - 1
  
  metaData <- unlist(strsplit(fname, split="\\."))
  
  if(length(metaData) == 8)
  {
    redundancy <- metaData[2]
    nGenes <- metaData[3]
    corrThresh <- metaData[5]
    candidate <- "None"
    filtering <- "None"
  }
  
  if(length(metaData) == 9)
  {
    if("MAD" == metaData[6])
    {
      redundancy <- metaData[2]
      nGenes <- metaData[3]
      corrThresh <- metaData[5]
      candidate <- "None"
      filtering <- "MAD"
    }
    else if("Quantile" == metaData[6])
    {
      redundancy <- metaData[2]
      nGenes <- metaData[3]
      corrThresh <- metaData[5]
      candidate <- "None"
      filtering <- "Quantile"
    }
    else
    {
      redundancy <- metaData[3]
      nGenes <- metaData[4]
      corrThresh <- metaData[6]
      candidate <- metaData[1]
      filtering <- "None"
    }
    
  }
  if(length(metaData) == 10)
  {
    redundancy <- metaData[3]
    nGenes <- metaData[4]
    corrThresh <- metaData[6]
    candidate <- metaData[1]
    filtering <- "MAD"

    if(metaData[7] == "Quantile")
    {
      filtering <- "Quantile"
    }
  }
  
  
  tmp$Redundancy <- redundancy
  tmp$Filtering <- filtering
  tmp$NumberMeasuredGenes <- nGenes
  tmp$TrueNumberMeasuredGenes <- trueNGenes
  tmp$CorrelationThreshold <- paste("0.",corrThresh,sep="")
  tmp$Candidates <- candidate
  tmp$GeneName <- rownames(tmp)
  return(tmp)
})

imputationResults.combined <- do.call(rbind, imputationResults)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#COMBINE COVERAGE AND IMPUTATION RESULTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

summaryTable <- getSummary(coverageResults, imputationResults.combined)

write.table(summaryTable, "Data/SummaryTable.Quantile.tsv",sep="\t",row.names=F,quote=F)




