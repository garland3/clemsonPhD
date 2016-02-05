__author__ = 'Anthony G'

from vtk import *
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
from vtk import vtkXMLUnstructuredGridWriter

from Node import Node
import os


class NodeHelper:
    pointDataVTK = []
    cellsVTK = []
    cellTypesVtk = []
    cellLocationsArray = []

    def calculateNodeSensitivity(self, listofElements, listOfNodes):
        # loop over the elements
        # get a list of the nodes for each element
        # add the element sensitivity to the node sensitivity and increment the number of elements contributing to a node each time we add to it

        # loop over the nodes and divide the sensitivity by the number of elements contributint to it.

        for e in listofElements:
            elementNodeList = e.nodesIndex  # list of nodes
            # print 'node list ', elementNodeList
            for i in elementNodeList:
                currentNum = int(
                    i) - 1  # minus 1 becasue the node numbers start at 1, but the list starts at 0, so subtract one to prevent overflow
                # print 'cur num' ,currentNum
                currentNode = listOfNodes[currentNum]
                currentNode.sensitivity += e.sensitivity
                currentNode.numberOfElements += 1

        c = 0
        print 'node sensitivty'
        for node in listOfNodes:
            c +=1
            numElements = node.numberOfElements
            node.sensitivity = node.sensitivity / float(numElements)
            print "{0} node elem {1}".format(node.sensitivity , numElements)

        return listOfNodes


    def calculateNodeDensities(self, listofElements, listOfNodes):
        # loop over the elements
        # get a list of the nodes for each element
        # add the element sensitivity to the node density and increment the number of elements contributing to a node each time we add to it

        # loop over the nodes and divide the density by the number of elements contributint to it.

        for e in listofElements:
            # print e.densityDesignVar
            elementNodeList = e.nodesIndex  # list of nodes
            # print 'node list ', elementNodeList
            for i in elementNodeList:
                currentNum = int(
                    i) - 1  # minus 1 becasue the node numbers start at 1, but the list starts at 0, so subtract one to prevent overflow
                # print 'cur num' ,currentNum
                currentNode = listOfNodes[currentNum]
                currentNode.density += e.densityDesignVar
                currentNode.numberOfElements += 1

        for node in listOfNodes:
            numElements = node.numberOfElements
            node.density = node.density / float(numElements)
            # print node.sensitivity

        return listOfNodes

    def readVTUfile(self,path, file_name):
        nodeList = []
        f = os.path.join(path, file_name)

        # Read the source file.
        reader = vtkXMLUnstructuredGridReader()
        reader.SetFileName(f)
        reader.Update()  # Needed because of GetScalarRange

        # Get the coordinates of nodes in the mesh
        nodes_vtk_array = reader.GetOutput().GetPoints().GetData()

        self.pointDataVTK = reader.GetOutput().GetPoints()
        self.cellsVTK = reader.GetOutput().GetCells()  # http://www.vtk.org/doc/release/4.2/html/classvtkUnstructuredGrid.html#z784_7
        self.cellTypesVtk = reader.GetOutput().GetCellTypesArray()  # http://www.vtk.org/doc/release/4.2/html/classvtkUnstructuredGrid.html#a4
        self.cellLocationsArray = reader.GetOutput().GetCellLocationsArray()  # http://www.vtk.org/doc/release/4.2/html/classvtkUnstructuredGrid.html#a5


        # Get the coordinates of the nodes and their displacements
        nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array)
        xlist, ylist, zlist = nodes_nummpy_array[:, 0], nodes_nummpy_array[:, 1], nodes_nummpy_array[:, 2]

        # The "Displacement" field is the first and only scalar in my vtk file
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
            newNode.setLocation(xlist[count], ylist[count], zlist[count])
            newNode.setDisplacement(disp[count, 0], disp[count, 1], disp[count, 2])
            nodeList.append(newNode)
            count += 1

        return nodeList


    def writeSensitivityVTKFile(self, rootFileDirectory,filename,listofElements, listOfNodes):
        f = os.path.join(rootFileDirectory,filename)
        da = vtkDoubleArray()
        da.SetName('sensitivity')
        da.SetNumberOfComponents(1)
        da.SetNumberOfTuples(len(listofElements))
        count = 1
        # with open('sensitivity.csv', "wt") as fout:
        # fout.write('x coord, y coord, z coord, scalar\n')
        for node in listOfNodes:
            da.InsertTuple1(count, node.sensitivity)
            count += 1
            # fout.write('{0}, {1} ,{2},{3}\n'.format(node.x,node.y,node.z,node.sensitivity))

        ugd = vtkUnstructuredGrid()
        ugd.SetPoints(self.pointDataVTK)

        # http://www.vtk.org/doc/release/4.2/html/classvtkUnstructuredGrid.html#z784_6
        ugd.SetCells(self.cellTypesVtk, self.cellLocationsArray, self.cellsVTK)

        upd = ugd.GetPointData()  # type(spd) is vtkPointData
        upd.SetScalars(da)

        # # http://www.programcreek.com/python/example/65341/vtk.vtkXMLUnstructuredGridWriter
        # # http://www.vtk.org/Wiki/VTK/Examples/Python/DataManipulation/Cube.py

        writer = vtkXMLUnstructuredGridWriter()
        writer.SetFileName(f)
        writer.SetDataModeToAscii()
        # # apparently some api stuff was changed
        # # http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Replacement_of_SetInput
        writer.SetInputData(ugd)
        writer.Update()
        writer.Write()


    def writeDensityVTKFile(self, folderpath, filename, listofElements, listOfNodes):



        densityData = vtkDoubleArray()
        densityData.SetName('density')
        densityData.SetNumberOfComponents(1)
        densityData.SetNumberOfTuples(len(listofElements))
        count = 1
        # with open('sensitivity.csv', "wt") as fout:
        # fout.write('x coord, y coord, z coord, scalar\n')
        for node in listOfNodes:
            densityData.InsertTuple1(count, node.density)
            count += 1
            # fout.write('{0}, {1} ,{2},{3}\n'.format(node.x,node.y,node.z,node.sensitivity))

        ugd = vtkUnstructuredGrid()
        ugd.SetPoints(self.pointDataVTK)

        # http://www.vtk.org/doc/release/4.2/html/classvtkUnstructuredGrid.html#z784_6
        ugd.SetCells(self.cellTypesVtk, self.cellLocationsArray, self.cellsVTK)

        upd = ugd.GetPointData()  # type(spd) is vtkPointData
        upd.SetScalars(densityData)


        # # http://www.programcreek.com/python/example/65341/vtk.vtkXMLUnstructuredGridWriter
        # # http://www.vtk.org/Wiki/VTK/Examples/Python/DataManipulation/Cube.py

        writer = vtkXMLUnstructuredGridWriter()
        f = os.path.join(folderpath,filename)
        writer.SetFileName(f) #"general.vtu")
        writer.SetDataModeToAscii()
        # # apparently some api stuff was changed
        # # http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Replacement_of_SetInput
        writer.SetInputData(ugd)

        writer.Update()
        writer.Write()


    def writeVTUFile(self, listofElements, listOfNodes):
        da = vtkDoubleArray()
        da.SetName('density')
        da.SetNumberOfComponents(1)
        da.SetNumberOfTuples(len(listofElements))
        count = 1
        # with open('sensitivity.csv', "wt") as fout:
        # fout.write('x coord, y coord, z coord, scalar\n')
        for node in listOfNodes:
            da.InsertTuple1(count, node.density)
            count += 1
            # fout.write('{0}, {1} ,{2},{3}\n'.format(node.x,node.y,node.z,node.sensitivity))

        ugd = vtkUnstructuredGrid()
        ugd.SetPoints(self.pointDataVTK)

        # http://www.vtk.org/doc/release/4.2/html/classvtkUnstructuredGrid.html#z784_6
        ugd.SetCells(self.cellTypesVtk, self.cellLocationsArray, self.cellsVTK)

        upd = ugd.GetPointData()  # type(spd) is vtkPointData
        upd.SetScalars(da)

        # # http://www.programcreek.com/python/example/65341/vtk.vtkXMLUnstructuredGridWriter
        # # http://www.vtk.org/Wiki/VTK/Examples/Python/DataManipulation/Cube.py

        writer = vtkXMLUnstructuredGridWriter()
        writer.SetFileName("density.vtu")
        writer.SetDataModeToAscii()
        # # apparently some api stuff was changed
        # # http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Replacement_of_SetInput
        writer.SetInputData(ugd)
        writer.Update()
        writer.Write()