import os


class Element:
	def __init__(self,number, body, nodes):
		self.nodeNumber = number
		self.bodyNumber = body
		self.nodes = nodes
		self.numberNodes = len(nodes)
		self.massMatrix = 0
		self.stiffnessMatrix = 0
		

		
		
# The main executable for this program.
def main():
	count=1
	elementsList = []
	rootFileDirectory = """C:\Users\Anthony G\Git\clemsonPhD\Experimental\ElmerFEA\\1TestProblemSIMP\sw\Mesh_1"""
	filenameInput = 'original_mesh.elements'
	filenameOUTput = 'mesh.elements'
	# ----------------------------
	# Read the existing elements file
	# ----------------------------
	f = os.path.join(rootFileDirectory,filenameInput)
	with open(f, "rt") as fin:
		for line in fin:
			values = [int(s) for s in line.split() if s.isdigit()]
			#elementNumber = re.search('/^[^\d]*(\d+)/', line)
			#print values
			l = len(values)
			g = Element(values[0], values[1], values[2:l])
			#print g.nodes
			elementsList.append(g)           
			
			count += 1

	# ----------------------------
	# Change all the elements to body 10
	# ----------------------------
	for element in elementsList:
		element.bodyNumber= 5
	

	# ----------------------------
	# Export the new elements file. 
	# ----------------------------
	ff = os.path.join(rootFileDirectory,filenameOUTput)
	with open(ff, "wt") as fout:
		for element in elementsList:
			#print element.nodes
			nodeString =  " ".join(map(str, element.nodes))
			#print nodeString
			#print element.nodeNumber
			text = "{0} {1} {2}\n".format(element.nodeNumber, element.bodyNumber, nodeString)
			#print text
			fout.write(text)
		
			
			
			
if __name__ == "__main__":
    main()
			