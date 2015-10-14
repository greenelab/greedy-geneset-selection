#!/bin/bash


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Global Variables
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#File to hold the error log
err="GGS.error.log"

#binary covariance matrix threshold
BinTHRESH="0.70 0.65 0.60"

#number of genes to select
NGENES="400 10 15 20 25 30 35 40 45 50 100 150 200 250 300 350"

#number of correlated genes selected in order to consider covered
NFOLD="3 2 1"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Identify the DM and Predictable sets
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for thresh in `echo $BinTHRESH`
do

	for fold in `echo $NFOLD`
	do

		for genes in `echo $NGENES`
		do
			echo Using BrCa.RNASeq.cor.$thresh.txt to Analyze $genes Genes at $fold fold coverage...
			../GGS.py BrCa.RNAseq.filtered.cor.$thresh.txt --fold=$fold --nGenes=$genes > GGS.Parameter.Sweep.Results/BrCa.TCGA.filtered.$fold.$genes.$thresh.results.txt 2>> GGS.error.log
			grep ^Final GGS.Parameter.Sweep.Results/BrCa.TCGA.filtered.$fold.$genes.$thresh.results.txt | cut -d " " -f 3- | sed -e 's/[ \t]*$//' | sed -e 's/\s/","/g' | sed -e 's/^/"/' | sed -e 's/$/"/' > GGS.Parameter.Sweep.Results/genesets/BrCa.TCGA.filtered.$fold.$genes.$thresh.measured.txt 2>> GGS.error.log
			grep ^"Covered Genes:" GGS.Parameter.Sweep.Results/BrCa.TCGA.filtered.$fold.$genes.$thresh.results.txt | cut -d " " -f 3- | sed -e 's/[ \t]*$//' | sed -e 's/\s/","/g' | sed -e 's/^/"/' | sed -e 's/$/"/' > GGS.Parameter.Sweep.Results/genesets/BrCa.TCGA.filtered.$fold.$genes.$thresh.covered.txt 2>> GGS.error.log


		done
	done

done
