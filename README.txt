Greedy Geneset Selection
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
Description:
This repository houses methods used to identify 'tag' genes in 
expression data such that a reduced number of genes is directly 
assayed but the expression values for a maximized number of 
genes can be predicted. We created a greedy approach to  
identify genes to directly measure which maximize the number of
genes whose expression we can predict. 


The greedy geneset selection (GGS) algorithm requires three arguments:
     1. A binary gene by gene symmetrical matrix in which 1 at position
   	i,j indicates a strong correlation between gene_i and gene_j
     2. Number of genes to directly measure
     3. Fold redundancy. Redundancy is the number of directly measured
   	genes that have to be correlated with an unmeasured gene in 
   	order for the unmeasured gene to be predictable. 

GGS can also be used in the context of existing candidate gene 
lists. Candidate genes are automatically added to the directly
measured geneset. 

For an explanation of the command line arguments use
python GGS.py -h 


GGS has been tested using Python 2.7 on Ubuntu 14.04
Used python libraries:
  docopt 0.6.1
  numpy 1.7.1
  itertools
  collections
  copy



Also included in this repository is the workflow used for analyses
and figures for the GGS manuscript (excluding those using non-public 
data from Konecny et al.) The bash script GGS.sh runs analyses 
sequentially and produces plots. 

Analyses were performed using R version 3.1.2
Used R libraries:
  ggplot2 1.0.0
  reshape 0.8.5
  survcomp 1.16.0
  survival 2.37
  randomForest 4.6
  outliers 0.14
  reshape2 1.4.1
  boot 1.3
  doMC 1.3.3
  ff 2.2
  curatedOvarianData 1.3.4
  igraph 0.7.1
  doppelgangR (https://github.com/lwaldron/doppelgangR)
    For convenience, scripts are included in the GGS.sh bash file to help install doppelgangR but are commented out. It is recommended that you install the latest version from the github repository but the included scripts worked as of 6/18/15

