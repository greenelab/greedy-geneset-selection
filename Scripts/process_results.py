#!/usr/bin/env python
"""
Usage:
	process_results.py <relativeDirectory> <matchingString> 
	process_results.py -h | --help


Options:
	-h --help 			  Show this screen.


"""

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LOAD PACKAGES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import glob
from docopt import docopt
import os
import errno

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#FUNCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#http://stackoverflow.com/questions/10840533/most-pythonic-way-to-delete-a-file-which-may-not-exist
def silentRemove(fName):
    try:
      os.remove(fName)
    except OSError as e:
      if e.errno != errno.ENOENT:
        raise
 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#GLOBAL VARIABLES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
VERSION="0.0.1"




summary_info = []

#Get which files to process


arguments = docopt(__doc__, version=VERSION)
matchString = arguments['<matchingString>']
relativeDirectory = arguments['<relativeDirectory>']

os.chdir(relativeDirectory)
oFName = "process_results." + matchString
silentRemove(oFName)
#Loop through all of the GGS result files
#for this_file in glob.glob("*65.results.txt"):
for this_file in glob.glob("*" + matchString):
	print(this_file)
	
	fakeNumMeasuredGenes = this_file.split(".")[2]
	measured_genes = None
	covered_genes = None
	fold = 0
	for line in open(this_file):


		#Get the selected Gene List
		if line.startswith('Final Genes: '):
			measured_genes = set(line.strip().split(" ")[2:])  #Start with index 2 to skip the 'Final Genes: '
			print("\t...Found Measured Genes")
			

		#Get the covered Gene List
		elif line.startswith('Covered Genes: '):
			covered_genes = set(line.strip().split(" ")[2:])  #Start with index 2 to skip the 'Covered Genes: '
			print("\t...Found Covered Genes")

		#Get the fold coverage
		elif line.startswith('Fold Coverage: '):
			fold = int(line.strip().split(" ")[2])

	if(measured_genes != None and covered_genes != None):
		additional_genes_covered = covered_genes.difference(measured_genes)

		summary_info.append([len(measured_genes), fold, len(covered_genes), len(additional_genes_covered), fakeNumMeasuredGenes])



print("======================================")
print("nMeasured\tFold\tScore\tnCovered\tFakeNMeasured")
for item in summary_info:
	for j in item:
		print(str(j) + "\t"),
	print("")

#oFile = open("process_results.65.txt", 'w')
oFile = open("process_results."+matchString, 'w')
oFile.write("nMeasured\tFold\tScore\tnCovered\tFakeNMeasured\tBlank\n")
for item in summary_info:
	for j in item:
		oFile.write(str(j) + "\t"),
	oFile.write("\n")