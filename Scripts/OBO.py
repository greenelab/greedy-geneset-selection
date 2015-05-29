import requests

"""Processes the OBO file containing the generic Gene Ontology (GO) Slim
list of pathways for homo sapiens. For biological process pathway in the
GO-slim we use tribe.greenelab.com to get the genes that belong to the 
pathway. The R script Eligible.R then uses the pathway and membership 
genes to perform enrichment. 

"""


# Define where Tribe is located.
TRIBE_URL = "http://tribe.greenelab.com"

#File which holds the OBO format of the Go-slim human pathway list
fname = "Data/goslim_generic.obo.txt"

#This class will be used to hold the key pathway info as the OBO file is parsed. 
class OBO:
	GO_namespace = None
	GO_name = None
	GO_ID = None
	GO_description = None

	def __init__(self, id):
		self.GO_ID = id



#Number of processes in the OBO file
numberProcesses = 0

#Number of processes belonging to Biological Processes
BP_counter = 0

#A list of the actual biological process OBO objects
BP_classes = []

#open an OBO file and iterate through it line by line
f = open(fname)

#Read the first line and process it
line = f.readline()
	
#Iterate through the remaining lines in the OBO file
while line:
	if line.strip() == "[Term]":
		numberProcesses = numberProcesses + 1
		line = f.readline()
		if(line.strip()[0:3] == "id:"):
			thisID = line.strip()[4:]
			thisOBO = OBO(thisID)
			line = f.readline()
			if(line.strip()[0:5] == "name:"):
				thisOBO.GO_name = line.strip()[6:]
				line = f.readline()
				if(line.strip()[0:10] == "namespace:"):
					thisOBO.GO_namespace = line.strip()[11:]
					line = f.readline()
					if(line.strip()[0:4] == "def:"):
						thisOBO.GO_description = line.strip()[5:]
						if(thisOBO.GO_namespace == "biological_process"):
							BP_classes.append(thisOBO)
							BP_counter = BP_counter + 1

	line = f.readline()


print("There are a total number of " + str(numberProcesses) + " terms in the generic GO Slim")
print(str(BP_counter) + " terms belong to Biological processes (" + str(100 * float(BP_counter) / numberProcesses) + "%)")

print("The following are the BP specific GO IDs which are in the generic GO Slim")
with open('Data/OBO.genelist.csv', 'w') as fp:
	for c in BP_classes:
		print(str(c.GO_ID) + " \t " + str(c.GO_name))

		#build the tribe request
		parameters = {'query': 'GO-BP-' + c.GO_ID[3:],
				  'show_tip': 'true',
				  'organism__scientific_name': 'Homo sapiens',
				  'xrdb': 'Symbol'}

		#query the tribe server
		r = requests.get(TRIBE_URL + '/api/v1/geneset/', params=parameters)
		result = r.json()

		#process the json result
		if len(result['objects']) > 0:
			tmp = result['objects'][0]
			tip = tmp['tip']
			genes = tip['genes']

			#create a ; seperated file with the process ID,
			#process name, and comma seperated gene list
			fp.write(c.GO_ID + ";" + c.GO_name + ";")
			thisString = str(genes[0])
			if len(genes) > 1:
				for x in genes[1:]:
					thisString = thisString + ", " + str(x)
			thisString = thisString + "\n"

			#write the string to the output file
			fp.write(thisString)

	
