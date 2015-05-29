#!/bin/bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Global variables
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#File to hold the error log
err="GGS.error.log"

#binary covariance matrix threshold
BinTHRESH="0.6 0.65 0.7"

#number of genes to select
NGENES="10 15 20 25 30 35 40 45 50 100 150 200 250 300 350 400"

#number of correlated genes selected in order to consider covered
NFOLD="3 2 1"

IMPUTATIONLIST=""
CANDIDATEIMPUTATIONLIST=""

YoshGeneList=`cat Data/Candidate.Genes/Yoshihara2012.sig.mapped.csv`
TCGAGeneList=`cat Data/Candidate.Genes/TCGA2011.sig.mapped.csv`


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Inclusion/Exclusion of samples
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
echo "Creating Inclusion/Exclusion list of samples..."
Rscript Scripts/Inclusion_doppelgangR_v2.R

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Data processing
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Create the TCGA Affymetrix expression correlation matrices and plots
echo "Creating Binary Correlation Matrices..."
Rscript Scripts/Process_Data_args.R $BinTHRESH 2>> R.error.log

#Define genes eligible to either be directly meaured or predicted
Rscript Scripts/Eligible.R

#Perform Enrichment on the eligible genes for |rP| = 0.6 and redundancy = 1


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Perform the Naive Gene selection using DM 
#size of 400 and 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for thresh in `echo $BinTHRESH`
do
	for fold in `echo $NFOLD`
	do
		echo Performing Naive gene selection using $thresh threshold and $fold redundancy...
		Rscript Scripts/Naive_args.R $thresh $fold 400 2>/dev/null > Data/Naive/Naive.$thresh.$fold.result.txt
	done	
done


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Parameter Sweep: Identify DM and Predictable
#sets using GGS accross the Correlation 
#Threshold, Redundancy, and number of DM genes
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for thresh in `echo $BinTHRESH`
do

	for fold in `echo $NFOLD`
	do

		for genes in `echo $NGENES`
		do
			
			echo Using TCGA.Quantile.cor.bin.$thresh.txt to Analyze $genes Genes at $fold fold coverage...
			./GGS.py Data/TCGA.Quantile.cor.bin.$thresh.txt --fold=$fold --nGenes=$genes > Data/Quantile.GGS.Parameter.Sweep.Results/TCGA.$fold.$genes.$thresh.results.txt 2>> GGS.error.log
			grep ^Final Data/Quantile.GGS.Parameter.Sweep.Results/TCGA.$fold.$genes.$thresh.results.txt | cut -d " " -f 3- | sed -e 's/[ \t]*$//' | sed -e 's/\s/","/g' | sed -e 's/^/"/' | sed -e 's/$/"/' > Data/Quantile.GGS.Parameter.Sweep.Results/genesets/TCGA.$fold.$genes.$thresh.measured.txt 2>> GGS.error.log
			grep ^"Covered Genes:" Data/Quantile.GGS.Parameter.Sweep.Results/TCGA.$fold.$genes.$thresh.results.txt | cut -d " " -f 3- | sed -e 's/[ \t]*$//' | sed -e 's/\s/","/g' | sed -e 's/^/"/' | sed -e 's/$/"/' > Data/Quantile.GGS.Parameter.Sweep.Results/genesets/TCGA.$fold.$genes.$thresh.covered.txt 2>> GGS.error.log

			IMPUTATIONLIST="$IMPUTATIONLIST TCGA.$fold.$genes.$thresh"

			
 			#None of the candidate gene lists are smaller than 100 genes, so only attempt to use them if DM size > 100
			if [ "$genes" -ge "100" ]
			then
				#GGS on Quantile TCGA with Yoshihara Gene List
				echo Using Yoshihara Candidate Genes and TCGA.Quantile.cor.bin.$thresh.txt to Analyze $genes Genes at $fold fold coverage...
				./GGS.py Data/TCGA.Quantile.candidates.cor.bin.$thresh.txt --fold=$fold --nGenes=$genes --candidates=$YoshGeneList > Data/Quantile.GGS.Parameter.Sweep.Results/Yoshihara.Candidate.Genes/TCGA.$fold.$genes.$thresh.results.txt 2>> GGS.error.log
				grep ^Final Data/Quantile.GGS.Parameter.Sweep.Results/Yoshihara.Candidate.Genes/TCGA.$fold.$genes.$thresh.results.txt | cut -d " " -f 3- | sed -e 's/[ \t]*$//' | sed -e 's/\s/","/g' | sed -e 's/^/"/' | sed -e 's/$/"/' > Data/Quantile.GGS.Parameter.Sweep.Results/Yoshihara.Candidate.Genes/genesets/TCGA.$fold.$genes.$thresh.measured.txt 2>> GGS.error.log
				grep ^"Covered Genes:" Data/Quantile.GGS.Parameter.Sweep.Results/Yoshihara.Candidate.Genes/TCGA.$fold.$genes.$thresh.results.txt | cut -d " " -f 3- | sed -e 's/[ \t]*$//' | sed -e 's/\s/","/g' | sed -e 's/^/"/' | sed -e 's/$/"/' > Data/Quantile.GGS.Parameter.Sweep.Results/Yoshihara.Candidate.Genes/genesets/TCGA.$fold.$genes.$thresh.covered.txt 2>> GGS.error.log

				
				#GGS on Quantile TCGA with TCGA Gene List
				echo Using TCGA Candidate Genes and TCGA.Quantile.cor.bin.$thresh.txt to Analyze $genes Genes at $fold fold coverage...
				./GGS.py Data/TCGA.Quantile.candidates.cor.bin.$thresh.txt --fold=$fold --nGenes=$genes --candidates=$TCGAGeneList > Data/Quantile.GGS.Parameter.Sweep.Results/TCGA.Candidate.Genes/TCGA.$fold.$genes.$thresh.results.txt 2>> GGS.error.log
				grep ^Final Data/Quantile.GGS.Parameter.Sweep.Results/TCGA.Candidate.Genes/TCGA.$fold.$genes.$thresh.results.txt | cut -d " " -f 3- | sed -e 's/[ \t]*$//' | sed -e 's/\s/","/g' | sed -e 's/^/"/' | sed -e 's/$/"/' > Data/Quantile.GGS.Parameter.Sweep.Results/TCGA.Candidate.Genes/genesets/TCGA.$fold.$genes.$thresh.measured.txt 2>> GGS.error.log
				grep ^"Covered Genes:" Data/Quantile.GGS.Parameter.Sweep.Results/TCGA.Candidate.Genes/TCGA.$fold.$genes.$thresh.results.txt | cut -d " " -f 3- | sed -e 's/[ \t]*$//' | sed -e 's/\s/","/g' | sed -e 's/^/"/' | sed -e 's/$/"/' > Data/Quantile.GGS.Parameter.Sweep.Results/TCGA.Candidate.Genes/genesets/TCGA.$fold.$genes.$thresh.covered.txt 2>> GGS.error.log

				CANDIDATEIMPUTATIONLIST="$CANDIDATEIMPUTATIONLIST Yoshihara.TCGA.$fold.$genes.$thresh TCGA.TCGA.$fold.$genes.$thresh"

			fi
		done
	done

done


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Use process_results.py to combine the parameter 
#sweep results into summary text files.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for thresh in `echo $BinTHRESH`
do

	echo Summarizing results for the $thresh Correlation Threshold analyses in Quantile data...
	python Scripts/process_results.py Data/Quantile.GGS.Parameter.Sweep.Results/ $thresh.results.txt

	echo Summarizing Yoshihara Candidate results for the $thresh Correlation Threshold analyses in Quantile data...
	python Scripts/process_results.py Data/Quantile.GGS.Parameter.Sweep.Results/Yoshihara.Candidate.Genes/ $thresh.results.txt

	echo Summarizing TCGA Candidate results for the $thresh Correlation Threshold analyses in Quantile data...
	python Scripts/process_results.py Data/Quantile.GGS.Parameter.Sweep.Results/TCGA.Candidate.Genes/ $thresh.results.txt

done

#Plot the number of predictable genes and eligible genes from the parameter sweep
echo "Plotting the parameter sweep results..."
Rscript Scripts/Plot_Quantile_Results_args.R $BinTHRESH


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Parameter Sweep: Perform expression imputation
#accross the Correlation Threshold, Redundancy, 
#and number of DM genes
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#R script to perform the imputation experiment
echo "Performing Imputation Experiment..."
echo "====================================\n====================================\n"
echo $IMPUTATIONLIST
echo "====================================\n====================================\n"
echo "Quantile Filtered Imputation..."
Rscript Imputation/Scripts/Imputation_Quantile_args.R $IMPUTATIONLIST 2>> R.error.log

#R script to perform the imputation experiment using candidate genelists
echo "Performing Imputation Experiment with Candidate Gene Lists..."
echo "====================================\n====================================\n"
echo $CANDIDATEIMPUTATIONLIST
echo "====================================\n====================================\n"
echo "Quantile Filtered Candidate Imputation"
Rscript Imputation/Scripts/Imputation_Quantile_args.R $CANDIDATEIMPUTATIONLIST 2>> R.error.log

#R script to plot the results of the imputation experiment
echo "Plotting results of imputation experiment..."
Rscript Imputation/Scripts/Plot_Quantile_Figures.R 2>> R.error.log

#R script to create the result table containing number of predictable genes and average prediction accuracy in TCGA and validation datasets
echo "Summarizing Results..."
Rscript Scripts/Summary_Tables_Quantile.R

#Create network diagram
echo "Creating the GML graph of genes to be used in Cytoskape..."
Rscript Scripts/GeneGraph.R


