#!/usr/bin/env python
"""
Usage:
	GGS.py <CorrelationMatrix> [--fold=<FOLD>] [--nGenes=<NG>] [--candidates=<CG>] [--permutations=<PERM>] [--seed=<SEED>] [--cpus=<NPROCESSES>]
	GGS.py -h | --help
	GGS.py --version

Options:
	-h --help 			  Show this screen.
	--version 			  Show the version.
	--fold=<FOLD>			  The minimum number of genes in the solution that must correlate with a specific gene in order for it to be considered covered.
        --nGenes=<NG> 			  Number of genes in solution excluding candidate genes
	--candidates=<CG>		  Candidate Gene Set
	--permutations=<PERM>		The number of random gene sets to draw in order to estimate a p-value for the selected solution
	--seed=<SEED>			Random seed; needed to ensure reproducibility of permutation testing
	--cpus=<NPROCESSES>		  Number of independent processes to use for fitness evaluation
	--log=<FILE>		  File to which EGS will output evolution statistics

"""


# The input for the CorrelationMatrix command-line argument is a '.txt' file, which is written by the 
# 'Process_Data.R' script, and contains a tab-delimited binary matrix.
# The "candidate genes" argument is a comma-separated list of genes.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#GLOBAL VARIABLES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
VERSION = '0.0.1-alpha'
VERBOSE = True

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#LOAD PACKAGES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import math					#This provides the absoulte value function
import sys					#This library allows the output stream to be flushed
from docopt import docopt			#This library provides the command line interface
import random as rand                           #This library provides the random integer functions for mutation and crossover
import numpy as np                              #This Library provides the binary array
from itertools import izip as zip, count	#This library provides the izip which is used to dynamically build the sets used to calculate coverage
from time import time				#Used to calculate the runtime of the evolution vs runtime of data load
from collections import Counter                 #Used to count the frequency of a gene in the coverage list. If freq > k, gene is covered.
import copy

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#GLOBAL VARIABLES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

data = None				#Holds the correlation matrix
data_const = None		#An unchanged version of data
candidates = None		#Holds the indeces for the candidate genes supplied
arguments = None		#Holds the command line arguments.
global FOLD 
FOLD= None
permutations = 0	
global SEED 
SEED = 0

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#FUNCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def eval_func(gene_idx_list, candidate_coverage=None):
	# eval_func takes a list of 'best' gene indices which is returned by the 'simpleGreedy' method

	global data
	data = data_const

	total_idxs_covered = [] # This will be a list of gene indexes which are covered by the 'best' genes 
	                        # passed in the parameter and candidate genes, if any
	if candidate_coverage:
	    total_idxs_covered.extend(candidate_coverage)

	for val in gene_idx_list:
	    t = genes_covered(int(val))
	    total_idxs_covered.extend(t)

	if FOLD == 1:
		score = len(set(total_idxs_covered))

	else:
		cnt = Counter(total_idxs_covered) # This counter holds a dictionary-style object where 
                                       # it counts how many times each gene (index) is covered
		genes_fold_covered = set([j for j in cnt if cnt[j] >= FOLD]) | set(gene_idx_list[:]) # Union of 2 sets, where the 2nd
                                                                                               # set is just genes covering themselves
		#TODO: Should we also add candidates to genes_fold_covered?
		score = len(genes_fold_covered)

	return score


def genes_covered(geneIdx):
        global data
	# This method is passed the index of a gene in the corr. matrix, and
	# it returns a list of gene indexes that are covered by the gene being passed
	corVals = data.data[geneIdx]
	covered = []

	covered = np.where(corVals == 1) # Find cells where value is 1
	covered = list(covered[0]) # Get the first column/row in numpy array (which are the gene indexes), 
                                   # and turn it into a list
	return(covered)


class corMat:

	def __init__(self, fName):
		self.data = []
		self.header = None
		self.load_cor_matrix(fName)


	def load_cor_matrix(self,fName):
        # Here, fName is the location of the file with the correlation matrix, which is the tab-separated text file  
        # written by the 'Process_Data.R' code.
            try:
                with open(fName) as input:
                    matrix_header = input.readline() # Read first line, which has the headers
                    self.header = matrix_header.replace('"', '').strip().split('\t')

                    for line in input: # Because readline() was just called, this will start from the second line
                        matrix_line = np.array(line.strip().split('\t')[1:], dtype=np.int8)
                        self.data.append(matrix_line)

            except IOError:
                print("Error: can\'t find the input gene correlation file or read the data")


def print_runtime_params():
	global arguments
	print("Binary Correlation Matrix: " + str(arguments['<CorrelationMatrix>']))
	print("Fold Coverage: " + str(arguments['--fold']))
	print("CPUs: " + str(arguments['--cpus']))
	print("Number of Genes: " + str(arguments['--nGenes']))
	print("Candidate Genes: " + str(arguments['--candidates']))
	print("-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- ")


def removeEdges(idx):
	global data

	for i in xrange(len(data.data)):
        # Loop to zero out all the values in the column that has the index 'idx', which was passed
        # (This will replace 1's with 0's and 0's with 0's)
		data.data[i][idx] = 0

        # This line places zeros in the whole row that has the index 'idx', which was passed
	data.data[idx] = np.zeros(len(data.data), dtype=np.int8)

	return 0


def simpleGreedy(matrix_size, soln_set_size, candidate_gene_idxs=None):

	total_gene_idxs = range(matrix_size)
	solution_set = set([])
	solution_list = [] # Solution list variable is made to preserve order of elements picked
	fold_counter = Counter() # This counter holds a dictionary-style object where it counts how many
                               # times each gene (index) is covered

	if candidate_gene_idxs:
	    for idx in candidate_gene_idxs:
	        fold_counter.update(genes_covered(idx)) # This is passed the list of gene indexes  
                                            # covered by candidates, returned by 'genes_covered()'
	        solution_set.add(idx)
	        solution_list.append(idx)

	already_covered = set([j for j in fold_counter if fold_counter[j] >= FOLD]) 

	for idx in already_covered:
	    removeEdges(idx)

	print("Covered by candidates: " + str(len(already_covered)))
	while (len(solution_set) < soln_set_size):
	    number_of_genes_covered = [(len(set(genes_covered(t))), t) for t in total_gene_idxs]
	    number_of_genes_covered.sort(reverse=True)
	    j = 0    
	    best = number_of_genes_covered[j][1]

	    while best in solution_set:
	        j += 1
	        best = number_of_genes_covered[j][1]

	    solution_set.add(best)
	    solution_list.append(best)
	    fold_counter.update(genes_covered(best))
	    now_covered = set([j for j in fold_counter if fold_counter[j] >= FOLD])
	    newly_covered = now_covered - already_covered
	    already_covered = now_covered

	    for idx in newly_covered:
	        removeEdges(idx)

	    #print("Number of Genes Selected: " + str(len(solution_set)) + "                          \r"),
	    #sys.stdout.flush()
	    this_gene = data.header[best]
	    print("Gene: " + str(this_gene) + "    " + str(len(newly_covered))  )
	    sys.stdout.flush()


	print("Number of Genes Selected: " + str(soln_set_size))
	sys.stdout.flush()
	result_dict = {}
	result_dict['solution_set'] = solution_set
	result_dict['solution_list'] = solution_list
	result_dict['fold_counter'] = fold_counter

	return result_dict


def main_run(corr_matrix_location, requested_num_genes, candidate_genes, perms=0):

	#Output command line options
	print_runtime_params()

	t0 = time()

        global data
	data = corMat(corr_matrix_location)

        global data_const
	data_const = copy.deepcopy(data)

	if VERBOSE:
		print("Loaded Binary Correlation Matrix")
		print(str(time() - t0) + " seconds")
		sys.stdout.flush()

	resulting_genes = [] # This will be a list of resulting gene indexes found by 'simpleGreedy'
	candidate_coverage = None

	if candidate_genes:
		provided_candidate_names = candidate_genes.strip().replace(" ", "").split(',') # Adding some code in case gene list
                                                                     # is formatted differently and includes some spaces
		global candidates
                candidates = []
                not_found = []
                for x in provided_candidate_names:
                    try: # Use in case candidate name is not found in corr_matrix
                        candidates.append(data.header.index(x))
                    except:
                        not_found.append(x)
                print('Genes not found in CorrelationMatrix: %s' % not_found) # Print list of genes not found in correlation matrix to console

	greedy_result_dict = simpleGreedy(len(data.data), requested_num_genes, candidate_gene_idxs=candidates)
	resulting_genes = greedy_result_dict['solution_list']
	fold_counter = greedy_result_dict['fold_counter']
	score = len(set([j for j in fold_counter if fold_counter[j] >= FOLD]))
	print("-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- ")

	if VERBOSE:
		print("Completed Greedy Gene Selection")
		print(str(time() - t0) + " seconds")
		sys.stdout.flush()

	message = "Final Genes: "
	for g in resulting_genes:
		message += str(data.header[g]) + " "
	print(message)
	
	covered_genes = [data.header[j] for j in fold_counter if fold_counter[j] >= FOLD]
	measured_genes = [data.header[j] for j in resulting_genes]
	message2 = "Covered Genes: "
	for g in covered_genes:
		message2 += str(g) + " "
	print(message2)
	print("Score: " + str(score))

	score_true = len(set(covered_genes).difference(set(measured_genes)))
	print("True Score: " + str(score_true))

	if perms > 0:
		print("Performing permutation test...")
		rand.seed(SEED)
		res_list = []
		for p in xrange(perms):
			gene_set = rand.sample(xrange(len(data.data)), requested_num_genes)
			res = get_set_coverage(gene_set)
			res_list.append(res)
			percent_complete = float((float((p + 1) / perms)) * 100)
			print(str(percent_complete) + " %                  \r"),
			sys.stdout.flush()

		pval = float(len([x for x in res_list if res_list[x] >= score_true]) / perms)
		print("\nPermuted Scores: " + str(res_list))
		print("P-value: " + str(pval))

	return(0)


def get_set_coverage(candidate_gene_idxs):

	solution_set = set([])
	solution_list = [] # Solution list variable is made to preserve order of elements picked
	fold_counter = Counter() # This counter holds a dictionary-style object where it counts how many
                               # times each gene (index) is covered

	if candidate_gene_idxs:
	    for idx in candidate_gene_idxs:
	        fold_counter.update(genes_covered(idx)) # This is passed the list of gene indexes  
                                            # covered by candidates, returned by 'genes_covered()'
	        solution_set.add(idx)
	        solution_list.append(idx)

	already_covered = set([j for j in fold_counter if fold_counter[j] >= FOLD]) 
	coverage = len(already_covered.difference(candidate_gene_idxs))
	return(coverage)

if __name__ == "__main__": #This helps to ensure Windows compatability, reads docopt arguments when code is run from command line

    arguments = docopt(__doc__, version=VERSION)

    if arguments['--fold'] is None:
        arguments['--fold'] = 1

    #global FOLD
    FOLD = int(arguments['--fold'])
    corr_matrix = arguments['<CorrelationMatrix>']
    requested_num_genes = int(arguments['--nGenes'])
    candidate_genes = arguments['--candidates']


    if arguments['--permutations'] is not None:
    	permutations = int(arguments['--permutations'])

    	if arguments['--seed'] is not None:
    		#global SEED
    		SEED = int(arguments['--seed'])

    main_run(corr_matrix, requested_num_genes, candidate_genes, permutations)
