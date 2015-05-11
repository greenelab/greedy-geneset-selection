#!/bin/bash

#binary covariance matrices
BinMATRIX="TCGA.cor.bin.0.6.txt TCGA.cor.bin.0.66.txt TCGA.cor.bin.0.7.txt"


#number of genes to select
NGENES=400


#number of permutations
PERMS=1000

#number of correlated genes selected in order to consider covered
NFOLD="3 2 1"


for fold in `echo $NFOLD`
do

	for binmat in `echo $BinMATRIX`
	do
		echo Analyzing $genes Genes at $fold fold coverage...
		./GGS.py $binmat --fold=$fold --nGenes=$NGENES --permutations=$PERMS > nullDist.$fold.$NGENES.$binmat.result.txt  &
		

	done
done