#````````````````````````````````````
#LIBRARIES
#````````````````````````````````````
library(reshape)
library(ggplot2)
library(grid)

#````````````````````````````````````
#FUNCTIONS
#````````````````````````````````````
MAD <- function(x)
{
  thisMedian <- median(x, na.rm=TRUE)
  result <- median(abs(x - thisMedian), na.rm=TRUE)
  return(result)
}


getImputationData <- function(fDir, fPatt, RNAseq=FALSE)
{
  files <- list.files(path=fDir, pattern=fPatt, full.names=T, recursive=FALSE)
  results <- lapply(files, function(x) {
    tmp <- read.table(x)
    fname <- tail(unlist(strsplit(x, split="/")), n=1)
    
    metaData <- unlist(strsplit(fname, split="\\."))
    if(length(metaData) == 9)
    {
      redundancy <- metaData[2]
      nGenes <- metaData[3]
      corrThresh <- metaData[5]
      candidate <- "None"
    }
    if(length(metaData) == 10)
    {
      redundancy <- metaData[3]
      nGenes <- metaData[4]
      corrThresh <- metaData[6]
      candidate <- metaData[1]
    }
    
    
    tmp$Redundancy <- redundancy
    tmp$NumberMeasuredGenes <- nGenes
    tmp$CorrelationThreshold <- paste("0.",corrThresh,sep="")
    tmp$Candidates <- candidate
    
    if(candidate == "Yoshihara" & as.numeric(nGenes) < 150)
    {
      tmp$NumberMeasuredGenes <- 121
    }
    else if(candidate == "TCGA" & as.numeric(nGenes) < 200)
    {
      tmp$NumberMeasuredGenes <- 183
    }
    return(tmp)
  })

  results.combined <- do.call(rbind, results)
  if(! RNAseq)
  {
    results.simple <- melt(results.combined[,c("TCGA.Spearman", "Tothill.Spearman", "Yoshihara.Spearman", "Mayo.Spearman","Bonome.Spearman", "Redundancy", "CorrelationThreshold", "NumberMeasuredGenes", "Candidates")], id=c("Redundancy", "CorrelationThreshold", "NumberMeasuredGenes", "Candidates"))  
  }
  else
  {
    results.simple <- melt(results.combined[,c("TCGA.RNAseq.ln.Spearman", "TCGA.RNAseq.ln.testing.Spearman", "Redundancy", "CorrelationThreshold", "NumberMeasuredGenes", "Candidates")], id=c("Redundancy", "CorrelationThreshold", "NumberMeasuredGenes", "Candidates"))
  }
  return(results.simple)
}


#````````````````````````````````````
#LOAD THE DATA
#````````````````````````````````````
#Get the  imputation results
dat <- getImputationData("Imputation/Data", "*.Quantile.imputation.results.txt$")
dat.RNA <- getImputationData("Imputation/Data", "*.Quantile.imputation.results.txt$", TRUE)

datasets <- unlist(lapply(levels(dat$variable), function(x){
  tmp <- strsplit(x, "\\.")
  return(unlist(tmp[[1]])[1])
}))
datasets.RNA <- gsub("\\.Spearman", "", levels(dat.RNA$variable))

dat$Plottable.Datasetname <- factor(dat$variable, labels=datasets)
dat$Plottable.DatasetnameII <- factor(dat$Plottable.Datasetname, levels=c("TCGA", "Mayo","Yoshihara", "Tothill", "Bonome"))
dat.RNA$Plottable.Datasetname <- factor(dat.RNA$variable, labels=datasets.RNA)
dat.RNA$Plottable.DatasetnameII <- factor(dat.RNA$Plottable.Datasetname, labels=c("paste('TCGA RNAseq All Samples')", "paste('TCGA RNAseq Testing Partition')"))

dat$Plottable.Threshold <- factor(dat$CorrelationThreshold, labels=c(bquote(abs(r[P]) >= "0.60"), bquote(abs(r[P]) >= "0.65"), bquote(abs(r[P]) >= "0.70"))) 
dat.RNA$Plottable.Threshold <- factor(dat.RNA$CorrelationThreshold, labels=c(bquote(abs(r[P]) >= "0.60"), bquote(abs(r[P]) >= "0.65"), bquote(abs(r[P]) >= "0.70"))) 

#````````````````````````````````````
#Plot the no candidates results
#````````````````````````````````````

# fig <- ggplot(dat[ dat$Candidates == "None",], aes(as.numeric(NumberMeasuredGenes), value, color=Redundancy, shape=Redundancy)) +
#   theme_bw() +
#   stat_summary(fun.data="mean_cl_boot", geom="pointrange", size=4, width=0.4) +
#   facet_grid(Plottable.DatasetnameII ~ Plottable.Threshold, scales="free_y", labeller=label_parsed) +  
#   xlab("Number of Measured Genes") +
#   ylab("Correlation of Predicted and Actual Expression") +
#   theme(axis.text.x=element_text(size=45), axis.text.y=element_text(size=45),
#         axis.title.x=element_text(size=50, vjust=-3, face="bold"), axis.title.y=element_text(size=50, vjust=3, face="bold"),
#         legend.title=element_text(size=50, face="bold"), legend.text = element_text(size=30, face="bold"),
#         legend.key.height=unit(5,"line"), legend.key.width=unit(5,"line"),
#         plot.margin = unit(c(1,6,1,6), "cm"),
#         panel.border = element_rect(colour="black", fill=NA, size=5))
# 
# fig.tmp <- fig + theme(panel.grid.major = element_line(colour = "grey", size = 2, linetype = 'solid'),
#                        panel.margin = unit(1, "lines"),
#                        strip.background = element_rect(color="black", fill="black", size=1, linetype="solid"),
#                        strip.text.x = element_text(color="white", face="bold", size=65),
#                        strip.text.y = element_text(color="white", face="bold", size=65))
# ggsave("Imputation/Figures/Quantile.No.Candidates.imputation.results.png", fig.tmp, height=51, width=51, limitsize=FALSE)
# ggsave("Imputation/Figures/Quantile.No.Candidates.imputation.results.pdf", fig.tmp, height=51, width=51, limitsize=FALSE)

dat$Redundancy <- factor(dat$Redundancy, labels=c("  1          ", "  2          ", "  3          "))
fig <- ggplot(dat[ dat$Candidates == "None",], aes(as.numeric(NumberMeasuredGenes), value, color=Redundancy, shape=Redundancy)) +
  theme_bw() +
  stat_summary(fun.data="mean_cl_boot", geom="pointrange", size=4, width=0.4) +
  facet_grid(Plottable.DatasetnameII ~ Plottable.Threshold, labeller=label_parsed) +  
  xlab("Number of Measured Genes") +
  ylab("Correlation of Predicted and Actual Expression") +
  theme(axis.text.x=element_text(size=45), axis.text.y=element_text(size=45),
        axis.title.x=element_text(size=50, vjust=-3, face="bold"), axis.title.y=element_text(size=50, vjust=3, face="bold"),
        legend.title=element_text(size=50, face="bold"), legend.text = element_text(size=50, face="bold"),
        plot.margin = unit(c(1,6,11,6), "cm"),
        panel.border = element_rect(colour="black", fill=NA, size=5))

fig.tmp <- fig + theme(panel.grid.major = element_line(colour = "grey", size = 2, linetype = 'solid'),
                       panel.margin = unit(1, "lines"),
                       strip.background = element_rect(color="black", fill="black", size=1, linetype="solid"),
                       strip.text.x = element_text(color="white", face="bold", size=65),
                       strip.text.y = element_text(color="white", face="bold", size=65))

fig.tmp2 <- fig.tmp + theme(legend.direction = "horizontal", legend.position=c(0.5, -0.08),
                            legend.key=element_rect(size=10, colour="white"), legend.key.size=unit(13.5, "lines"))
ggsave("Imputation/Figures/Quantile.No.Candidates.imputation.results.pdf", fig.tmp2, height=51, width=51, limitsize=FALSE)

#````````````````````````````````````
#Plot the TCGA candidates results
#````````````````````````````````````

fig <- ggplot(dat[ dat$Candidates == "TCGA",], aes(as.numeric(NumberMeasuredGenes), value, color=Redundancy, shape=Redundancy)) +
  theme_bw() + 
  stat_summary(fun.data="mean_cl_boot", geom="pointrange", size=3, width=0.4) +
  facet_grid(Plottable.DatasetnameII ~ Plottable.Threshold, scales="free_y", labeller=label_parsed) +  
  xlab("Number of Measured Genes") +
  ylab("Correlation of Predicted and Actual Expression") +
  theme(axis.text.x=element_text(size=45), axis.text.y=element_text(size=45),
        axis.title.x=element_text(size=50, vjust=-3, face="bold"), axis.title.y=element_text(size=50, vjust=3, face="bold"),
        legend.title=element_text(size=50, face="bold"), legend.text = element_text(size=45, face="bold"),
        legend.key.height=unit(5,"line"), legend.key.width=unit(5,"line"),
        plot.margin = unit(c(1,6,1,6), "cm"),
        panel.border = element_rect(colour="black", fill=NA, size=5),
        panel.grid.major = element_line(colour = 'grey', size = 2, linetype = 'solid'),
        panel.margin = unit(1, "lines"),
        strip.background = element_rect(color="black", fill="black", size=0.1, linetype="solid"),
        strip.text.x = element_text(color="white", face="bold", size=65),
        strip.text.y = element_text(color="white", face="bold", size=65, vjust=1))
ggsave("Imputation/Figures/Quantile.TCGA.Candidates.imputation.results.png", fig, height=51, width=51, limitsize=FALSE)




#````````````````````````````````````
#Plot the Yoshihara candidates results
#````````````````````````````````````

fig <- ggplot(dat[ dat$Candidates == "Yoshihara",], aes(as.numeric(NumberMeasuredGenes), value, color=Redundancy, shape=Redundancy)) +  
  theme_bw() +
  stat_summary(fun.data="mean_cl_boot", geom="pointrange", size=3, width=0.4) +
  facet_grid(Plottable.DatasetnameII ~ Plottable.Threshold, scales="free_y", labeller=label_parsed) +  
  xlab("Number of Directly Measured Genes") +
  ylab("Correlation of Predicted and Actual Expression") +
  theme(axis.text.x=element_text(size=45), axis.text.y=element_text(size=45),
        axis.title.x=element_text(size=50, vjust=-3, face="bold"), axis.title.y=element_text(size=50, vjust=3, face="bold"),
        legend.title=element_text(size=50, face="bold"), legend.text = element_text(size=45, face="bold"),
        legend.key.height=unit(5,"line"), legend.key.width=unit(5,"line"),
        plot.margin = unit(c(1,6,1,6), "cm"),
        panel.border = element_rect(colour="black", fill=NA, size=5),
        panel.grid.major = element_line(colour = 'grey', size = 2, linetype = 'solid'),
        panel.margin = unit(1, "lines"),
        strip.background = element_rect(color="black", fill="black", size=0.1, linetype="solid"),
        strip.text.x = element_text(color="white", face="bold", size=65),
        strip.text.y = element_text(color="white", face="bold", size=65, vjust=1))
ggsave("Imputation/Figures/Quantile.Yoshihara.Candidates.imputation.results.png", fig, height=51, width=51, limitsize=FALSE)



 
#````````````````````````````````````
#Plot the RNAseq results
#````````````````````````````````````
fig <- ggplot(dat.RNA[dat.RNA$Candidates == "None",], aes(as.numeric(NumberMeasuredGenes), value, color=Redundancy, shape=Redundancy)) +  
  theme_bw() +
  stat_summary(fun.data="mean_cl_boot", geom="pointrange", size=3, width=0.4) +
  facet_grid(Plottable.DatasetnameII ~ Plottable.Threshold, labeller=label_parsed) +  
  xlab("Number of Directly Measured Genes") +
  ylab("Correlation of Predicted and Actual Expression") +
  theme(axis.text.x=element_text(size=45), axis.text.y=element_text(size=45),
        axis.title.x=element_text(size=50, vjust=-3, face="bold"), axis.title.y=element_text(size=50, vjust=3, face="bold"),
        legend.title=element_text(size=50, face="bold"), legend.text = element_text(size=45, face="bold"),
        legend.key.height=unit(5,"line"), legend.key.width=unit(5,"line"),
        plot.margin = unit(c(1,6,1,6), "cm"),
        panel.border = element_rect(colour="black", fill=NA, size=5),
        panel.grid.major = element_line(colour = 'grey', size = 2, linetype = 'solid'),
        panel.margin = unit(1, "lines"),
        strip.background = element_rect(color="black", fill="black", size=0.1, linetype="solid"),
        strip.text.x = element_text(color="white", face="bold", size=65),
        strip.text.y = element_text(color="white", face="bold", size=65, vjust=1))
ggsave("Imputation/Figures/Quantile.RNAseq.imputation.results.png", fig, height=31, width=51, limitsize=FALSE)


