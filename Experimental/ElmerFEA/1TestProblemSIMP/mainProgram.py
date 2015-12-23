from vtk import *
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

class ElementHelper:	
	def __init__(self):
		self.example = 1
		
	def readMeshFile(self, filename):
		count=1
		elementsList = []
		# ----------------------------
		# Read the existing elements file
		# ----------------------------
		with open(filename, "rt") as fin:
			for line in fin:
				values = [int(s) for s in line.split() if s.isdigit()]
				#elementNumber = re.search('/^[^\d]*(\d+)/', line)
				#print values
				
				
				l = len(values)
				newElement = Element(values[0], values[1], values[2:l])
				#print g.nodes
				elementsList.append(newElement)           
				
				count += 1
		return elementsList
		
	def writeMeshFile(self, filename,elementsList):
		
		# ----------------------------
		# Export the new elements file. 
		# ----------------------------
		with open(filename, "wt") as fout:
			for element in elementsList:
				#print element.nodes
				nodeString =  " ".join(map(str, element.nodesIndex))
				#print nodeString
				#print element.nodeNumber
				text = "{0} {1} {2}\n".format(element.nodeNumber, element.bodyNumber, nodeString)
				#print text
				fout.write(text)
				
	def readLocalStiffnessFile(self,filename,elementsList):
		with open(filename, "rt") as fin:
			count =0
			countElement = 0 # the element indexes start at 0?? I think??
			for line in fin:
				# REad 17 lines
				values = [int(s) for s in line.split() if s.isdigit()]
		

class Element:
	def __init__(self,number, body, nodesindex):
		self.nodeNumber = number
		self.bodyNumber = body
		self.nodesIndex = nodesindex
		self.numberNodes = 0
		self.massMatrix = 0
		self.stiffnessMatrix = 0
	
class Node:
	# the node number
	number = 0
	
	# location in 3 coordinates
	x = 0
	y = 0
	z = 0
	
	# displacement in the 3 dof
	disp_x = 0
	disp_y = 0
	disp_z = 0 	
	
	def __init__(self,Nodenumber):
		self.number = Nodenumber
	
	def setLocation(self,x,y,z):
		self.x = x
		self.y = y
		self.z = z
	
	def setDisplacement(self,dx,dy,dz):
		self.disp_x = dx
		self.disp_y = dy
		self.disp_z = dz
	
class NodeHelper:
	def readVTUfile(self, file_name):	
		nodeList=[]
		 
		# Read the source file.
		reader = vtkXMLUnstructuredGridReader()
		reader.SetFileName(file_name)
		reader.Update() # Needed because of GetScalarRange

		# Get the coordinates of nodes in the mesh
		nodes_vtk_array= reader.GetOutput().GetPoints().GetData()

		#Get the coordinates of the nodes and their displacements
		nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array)
		xlist,ylist,zlist= nodes_nummpy_array[:,0] , nodes_nummpy_array[:,1] , nodes_nummpy_array[:,2]
		
		#The "Displacement" field is the first and only scalar in my vtk file
		displacement_vtk_array = reader.GetOutput().GetPointData().GetArray(0)
		disp = vtk_to_numpy(displacement_vtk_array)
		
		print "number of nodes is " + str(len(xlist))
		
		# -----------
		# make the new nodes
		# -------------
		
		count = 0
		for x in xlist:
			#print i
			newNode = Node(count)
			newNode.setLocation(xlist[count],ylist[count],zlist[count])
			newNode.setDisplacement(disp[count,0],disp[count,1],disp[count,2])
			nodeList.append(newNode)
			
		
		return nodeList
			
		
	
		
		
# The main executable for this program.
def main():
	helper = ElementHelper()
	elementsList = helper.readMeshFile('original_mesh.elements')
	
	# Read the vtk file and get the node locations and displacements
	file_name = "paraviewout0001.vtu"
	nHelper = NodeHelper()
	nodeList = nHelper.readVTUfile(file_name)
	
	nHelper.readLocalStiffnessFile('lstiffness.txt',elementsList)
	
	
	# ----------------------------
	# Change all the elements to body 10
	# ----------------------------
	#for element in elementsList:
	#	element.bodyNumber= 5
	
	# helper.writeMeshFile('mesh.elements',elementsList)
	
		
			
			
			
if __name__ == "__main__":
    main()
			