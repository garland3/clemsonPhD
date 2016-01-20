from vtk import *
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
from vtk import vtkXMLUnstructuredGridWriter

from ElementHelper import ElementHelper
from Node import Node
from Element import Element
from NodeHelper import NodeHelper

from subprocess import call

import os, sys, shutil

# calculate the element sensitivity
def calculateElementSensitivity(elementsList, nodeList):
    # loop over the elements
    # the local sensitivity is
    # s = mass*u*K*u^T
    for e in elementsList:
        e.mass = sum(e.massMatrix)
        localDisplacementVector = []

        # get the element displacement matrix
        # 1. Figure out which nodes are in this element
        # 2. Get the displacements of those nodes and make a vector
        elementNodeList = e.nodesIndex  # list of nodes
        # print 'node list ', elementNodeList
        for i in elementNodeList:
            currentNum = int(i) - 1  # minus 1 becasue the node numbers start at 1,
                #  but the list starts at 0, so subtract one to prevent overflow
                #	print 'cur num' ,currentNum
            currentNode = nodeList[currentNum]
            localDisplacementVector.append(currentNode.disp_x)
            localDisplacementVector.append(currentNode.disp_y)
            localDisplacementVector.append(currentNode.disp_z)

        localDispVector = np.asarray(localDisplacementVector)  # convert to array
        localStiffness = np.asarray(e.stiffnessMatrix)  # convert to np array
        DOF = 3
        sizes = e.numberNodes * DOF  # should be 12

        resizedLStiff = np.reshape(localStiffness, [sizes, sizes])

        tr = np.transpose(localDispVector)

        penalty = 3
        # -penal*x(ely,elx)^(penal-1)*Ue'*KE*Ue;
        x= e.densityDesignVar
        s = - penalty * x**(penalty-1) * np.dot(np.dot(localDispVector, resizedLStiff), tr) * e.mass
        print s
        #s = (np.dot(np.dot(localDispVector, resizedLStiff), tr) * e.mass)**3
        # s = e.densityDesignVar #penalize intermediate values

        e.sensitivity = s*1e10  # store the sensitivity

def callElmer(loopNumber):
        cmd = "./ElmerSolver case.sif"
        call(cmd)


# The main executable for this program.
def main():
    rootFileDirectory = """C:\Users\Anthony G\Git\clemsonPhD\Experimental\ElmerFEA\\1TestProblemSIMP\sw\Mesh_1"""
    helper = ElementHelper()
    elementsList = helper.readMeshFile(rootFileDirectory,'original_mesh.elements')

    # Read the vtk file and get the node locations and displacements
    file_name = "paraviewout0001.vtu"
    nHelper = NodeHelper()
    nodeList = nHelper.readVTUfile(rootFileDirectory,file_name)
    helper.readLocalStiffnessFile(rootFileDirectory,'lstiffness.txt', elementsList)
    helper.readLocalMassFile(rootFileDirectory,'lmass.txt', elementsList)
    calculateElementSensitivity(elementsList, nodeList)
    # these 2 functions are for visualization only
    filename = 'sensitivity.vtu'
    nodeList = nHelper.calculateNodeSensitivity(elementsList, nodeList)
    nHelper.writeSensitivityVTKFile(rootFileDirectory,filename, elementsList, nodeList)

    targetVolFraction = 0.5
    helper.updateElementDensities(elementsList, nodeList,targetVolFraction)
    helper.writeElementMeshFile(rootFileDirectory,'mesh.elements', elementsList)

    # these 2 functions are for visualization only
    # Calcualte the node sensitivies so we can plot them in paraview
    nodeList = nHelper.calculateNodeDensities(elementsList, nodeList)
    nHelper.writeDensityVTKFile(rootFileDirectory,'density.vtu',elementsList, nodeList)


    # nHelper.writeDensityVTKFile(elementsList, nodeList)


if __name__ == "__main__":
    main()