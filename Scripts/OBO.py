


class OBO:
	GO_namespace = None
	GO_name = None
	GO_ID = None
	GO_description = None

	def __init__(self, id):
		self.GO_ID = id


fname = "Data/goslim_generic.obo.txt"

numberProcesses = 0
BP_counter = 0
BP_classes = []

#open an OBO file and iterate through it line by line
f = open(fname)

line = f.readline()

if line.strip() == "[Term]":
	numberProcesses = numberProcesses + 1
	line = f.readline()
	
	


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
for c in BP_classes:
	print(str(c.GO_ID) + " \t " + str(c.GO_name))
	#print(str(c.GO_name))
