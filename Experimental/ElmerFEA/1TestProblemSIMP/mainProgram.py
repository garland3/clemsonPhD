from vtk import *
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import re
from vtk import vtkXMLUnstructuredGridWriter
# from evtk.hl import pointsToVTK 

class ElementHelper:	
	def __init__(self):
		self.example = 1
		
	# read the mesh.elememnt file and get the element number, type, and node number array
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
				
				# number, body, typeNumber, nodesindex
				newElement = Element(values[0], values[1], values[2],values[3:l])
				newElement.numberNodes = l-3 # record the number of nodes
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
			count = 0
			countElement = -1 # the element indexes start at 0?? I think??
			
			currentElement  = elementsList[0]
			newElement = re.compile('--')
			
			# loop over the lines of the file
			for line in fin:
				
				# if found a '--' then a new element
				m = newElement.search(line)
				if m:
					found = m.group(0)
					countElement = countElement+1
					# print countElement
					
					# update the current element
					currentElement  = elementsList[countElement]
					count = 0
					continue  # go to the next line. 
					
				# add the current value to the stiffness matrix 
				currentElement.stiffnessMatrix.append(float(line))
				count = count+1
				
	def readLocalMassFile(self,filename,elementsList):
		with open(filename, "rt") as fin:
			count = 0
			countElement = -1 # the element indexes start at 0?? I think??
			
			currentElement  = elementsList[0]
			newElement = re.compile('--')
			
			# loop over the lines of the file
			for line in fin:
				
				# if found a '--' then a new element
				m = newElement.search(line)
				if m:
					found = m.group(0)
					countElement = countElement+1
					# print countElement
					
					# update the current element
					currentElement  = elementsList[countElement]
					count = 0
					continue  # go to the next line. 
					
				# add the current value to the mass matrix 
				currentElement.massMatrix.append(float(line))
				count = count + 1
	
	
class Element:
	def __init__(self,number, body, typeNumber, nodesindex):
		self.nodeNumber = number
		self.bodyNumber = body
		self.elementTypeNumber = typeNumber
		self.nodesIndex = nodesindex # list with nodes in this element
		self.numberNodes = 0
		self.massMatrix = []
		self.stiffnessMatrix = []
		self.mass = 0
		self.densityDesignVar = 1
		self.sensitivity = 0;
	
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
	
	# sensitivity metrics
	sensitivity = 0
	numberOfElements = 0 # the number of elements contributing to this node 
	
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
	
	pointDataVTK = []
	cellDataVtk = []
	cellDataTypesVTK =  []
	
	def readVTUfile(self, file_name):	
		nodeList=[]
		 
		# Read the source file.
		reader = vtkXMLUnstructuredGridReader()
		reader.SetFileName(file_name)
		reader.Update() # Needed because of GetScalarRange

		# Get the coordinates of nodes in the mesh
		nodes_vtk_array= reader.GetOutput().GetPoints().GetData()
		
		self.pointDataVTK = reader.GetOutput().GetPoints()
		self.cellDataVtk = reader.GetOutput().GetCells()
		self.cellDataTypesVTK = reader.GetOutput().GetCells().GetData()

		#Get the coordinates of the nodes and their displacements
		nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array)
		xlist,ylist,zlist= nodes_nummpy_array[:,0] , nodes_nummpy_array[:,1] , nodes_nummpy_array[:,2]
		
		#The "Displacement" field is the first and only scalar in my vtk file
		displacement_vtk_array = reader.GetOutput().GetPointData().GetArray(0)
		disp = vtk_to_numpy(displacement_vtk_array)
		
		print "number of nodes is " + str(len(xlist))
		# print disp
		
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
			count+=1
			
		
		return nodeList
		
	def writeSensitivityVTKFile(self,listofElements, listOfNodes):
		# loop over the elements
		# get a list of the nodes for each element
		# add the element sensitivity to the node sensitivity and increment the number of elements contributing to a node each time we add to it
		
		# loop over the nodes and divide the sensitivity by the number of elements contributint to it. 
				
		for e in listofElements:			
			elementNodeList = e.nodesIndex # list of nodes
			#print 'node list ', elementNodeList 
			for i in elementNodeList:
				currentNum = int(i) -1 # minus 1 becasue the node numbers start at 1, but the list starts at 0, so subtract one to prevent overflow
				#	print 'cur num' ,currentNum
				currentNode = listOfNodes[currentNum]		
				currentNode.sensitivity += e.sensitivity
				currentNode.numberOfElements +=1
		
		
		# p = vtkPoints()
		# p.SetNumberOfPoints(len(listofElements))
		# p.SetDataTypeToDouble()
		
		da = vtkDoubleArray()
		da.SetName('sensitivity')
		da.SetNumberOfComponents(1)
		da.SetNumberOfTuples(len(listofElements))
		
		
		#for i,pos in enumerate(atoms.get_positions()):
		#	p.InsertPoint(i,pos[0],pos[1],pos[2])
		
		# listOfXPoints = []
		# listOfYPoints = []
		# listOfZPoints = []
		listOfSensitivities = []
		count = 1
		for node in listOfNodes:
			# numElements = node.numberOfElements
			# node.sensitivity = node.sensitivity/float(numElements)
			# print node.sensitivity
			
			#listOfXPoints.append(node.x)
			#listOfYPoints.append(node.y)
			#listOfZPoints.append(node.z)
			#listOfSensitivities.append(node.sensitivity)
			# p.InsertPoint(count,node.x,node.y,node.z)
			da.InsertTuple1(count,node.sensitivity)
			count+=1
		
		ugd = vtkUnstructuredGrid()
		ugd.SetPoints(self.pointDataVTK)
		print self.cellDataTypesVTK
		#ugd.SetCells(self.cellDataVtk)
		
		upd = ugd.GetPointData() # type(spd) is vtkPointData
		upd.SetScalars(da)
		
		# http://www.programcreek.com/python/example/65341/vtk.vtkXMLUnstructuredGridWriter
		# http://www.vtk.org/Wiki/VTK/Examples/Python/DataManipulation/Cube.py
		
		# print vtkXMLUnstructuredGridWriter()
		#w = vtkXMLUnstructuredGridWriter()
		
		writer = vtkXMLUnstructuredGridWriter()
		writer.SetFileName("sensitivity.vtu")
		
		# p.SetScalars(scalars)
		#writer.GetCompressor().SetCompressionLevel(0)
		writer.SetDataModeToAscii()
		
		# apparently some api stuff was changed 
		# http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Replacement_of_SetInput
		writer.SetInputData(ugd)
		writer.Update()
		writer.Write()
		
        #writer.SetInput(self.Mesh)
        #w.SetFileName('sensitivity')
        
		# if self.Mode == "binary":
         #   writer.SetDataModeToBinary()
        #elif self.Mode == "ascii":
        #w.SetDataModeToAscii()
        #w.Write()
			
		#vtk.pointsToVTK("./points", listOfXPoints, listOfYPoints, listOfZPoints, data = {"sensitivity" : listOfSensitivities})

			
def create_points(array):
    """Create vtkPoints from double array"""
    vtk_points = vtkPoints()
    double_array = vtkDoubleArray()
    double_array.SetVoidArray(array, len(array), 1)
    double_array.SetNumberOfComponents(3)
    vtk_points.SetData(double_array)
    return vtk_points

	
# calculate the element sensitivity
def calculateElementSensitivity(elementsList,nodeList):
	# loop over the elements
	# the local sensitivity is
	# s = mass*u*K*u^T
	for e in elementsList:
		e.mass = sum(e.massMatrix)
		# print e.massMatrix
		# print e.mass
		
		localDisplacementVector = []
		
		# get the element displacement matrix
		# 1. Figure out which nodes are in this element
		# 2. Get the displacements of those nodes and make a vector
		elementNodeList = e.nodesIndex # list of nodes
		#print 'node list ', elementNodeList 
		for i in elementNodeList:
			currentNum = int(i) -1 # minus 1 becasue the node numbers start at 1, but the list starts at 0, so subtract one to prevent overflow
		#	print 'cur num' ,currentNum
			currentNode = nodeList[currentNum]
			localDisplacementVector.append(currentNode.disp_x)
			localDisplacementVector.append(currentNode.disp_y)
			localDisplacementVector.append(currentNode.disp_z)
			
		localDispVector = np.asarray(localDisplacementVector) # convert to array
		# print localDisplacementVector
		
		localStiffness = np.asarray(e.stiffnessMatrix) # convert to np array
		
		DOF = 3
		sizes = e.numberNodes*DOF # should be 12
		# print sizes
		
		#resizedLStiff= np.reshape(localStiffness,(sizes, -1)) # convert to 12 by 12 array
		#resizedLStiff = localStiffness.resize((sizes,sizes))
		resizedLStiff=np.reshape(localStiffness,[12,12])
		
		# print localStiffness
		
		
		tr = np.transpose(localDispVector)
		#print localDispVector.shape
		#print tr.size
		#print localStiffness.shape
		#print resizedLStiff.shape
		#print e.mass.size
		
		#s = localDispVector*resizedLStiff*tr*e.mass
		s = np.dot(np.dot(localDispVector,resizedLStiff),tr)*e.mass
		e.sensitivity = s # store the sensitivity
		print s
			
		
	
		
# The main executable for this program.
def main():
	helper = ElementHelper()
	elementsList = helper.readMeshFile('original_mesh.elements')
	
	# Read the vtk file and get the node locations and displacements
	file_name = "paraviewout0001.vtu"
	nHelper = NodeHelper()
	nodeList = nHelper.readVTUfile(file_name)
	
	helper.readLocalStiffnessFile('lstiffness.txt',elementsList)
	
	helper.readLocalMassFile('lmass.txt',elementsList)
	
	calculateElementSensitivity(elementsList,nodeList)
	
	nHelper.writeSensitivityVTKFile(elementsList,nodeList)
	
	#for n in nodeList:
	#	disp = []
	#	disp.append(n.disp_x)
	#	disp.append(n.disp_y)
	#	disp.append(n.disp_z)
	#	print disp
	
	
	# ----------------------------
	# Change all the elements to body 10
	# ----------------------------
	#for element in elementsList:
	#	element.bodyNumber= 5
	
	# helper.writeMeshFile('mesh.elements',elementsList)
	
		
			
			
			
if __name__ == "__main__":
    main()
			