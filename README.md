#Greedy Geneset Selection

#### By [James Rudd](http://www.dartmouth.edu/~doherty/personnel.html), [Rene A. Zelaya](http://www.greenelab.com/lab-members), and [Casey Greene](http://www.greenelab.com/)

[![DOI](https://zenodo.org/badge/18768/greenelab/greedy-geneset-selection.svg)](https://zenodo.org/badge/latestdoi/18768/greenelab/greedy-geneset-selection)

## Description

This repository houses methods used to identify 'tag' genes in expression data such that a reduced number of genes is directly assayed but the expression values for a maximized number of genes can be predicted. We created a greedy approach to identify genes to directly measure which maximize the number of genes whose expression we can predict.

The greedy geneset selection (GGS) algorithm requires three arguments:

1. A binary gene by gene symmetrical matrix in which 1 at position 	i,j indicates a strong correlation between gene<sub>i</sub> and gene<sub>j</sub>

2. Number of genes to directly measure

3. Fold redundancy. Redundancy is the number of directly measured genes that have to be correlated with an unmeasured gene in order for the unmeasured gene to be predictable. 

GGS can also be used in the context of existing candidate gene 
lists. Candidate genes are automatically added to the directly
measured geneset. 

For an explanation of the command line arguments use
python GGS.py -h 






Also included in this repository is the workflow used for analyses and figures for the GGS manuscript (excluding those using non-public data from Konecny et al.) The bash script GGS.sh runs analyses sequentially and produces plots. Seperate and confirmatory analyses were performed using TCGA RNAseq data and R scripts for these are located in the BrCa subdirectory. 

##Dependencies

GGS has been tested using Python 2.7 on Ubuntu 14.04

Used python libraries:

-  docopt 0.6.1
-  numpy 1.7.1
-  itertools
-  collections
-  copy



Analyses  for the publication were performed using R version 3.1.2

R libraries used:

-  ggplot2 1.0.0
-  reshape 0.8.5
-  survcomp 1.16.0
-  survival 2.37
-  randomForest 4.6
-  outliers 0.14
-  reshape2 1.4.1
-  boot 1.3
-  doMC 1.3.3
-  ff 2.2
-  curatedOvarianData 1.3.4
-  igraph 0.7.1
- RTCGAToolbox 1.99.4
-  doppelgangR (https://github.com/lwaldron/doppelgangR  last tested using commit sha: f0e02bca770efab4671a8290cd2580f0ec5cc674)
    

For convenience, scripts are included in the GGS.sh bash file to help install doppelgangR but are commented out. It is recommended that you install the latest version from the github repository but the included scripts worked as of 6/18/15
