__author__ = 'Anthony G'

from Element import Element
import re
import math
import os

from vtk import *
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
from vtk import vtkXMLUnstructuredGridWriter


class ElementHelper:
    def __init__(self):
        self.example = 1

    # read the mesh.elememnt file and get the element number, type, and node number array
    def readMeshFile(self, path, filename):
        count = 1
        elementsList = []
        # ----------------------------
        # Read the existing elements file
        # ----------------------------
        f = os.path.join(path, filename)
        print f
        with open(f, "rt") as fin:
            for line in fin:
                values = [int(s) for s in line.split() if s.isdigit()]
                # elementNumber = re.search('/^[^\d]*(\d+)/', line)
                #print values


                l = len(values)

                # number, body, typeNumber, nodesindex
                newElement = Element(values[0], values[1], values[1]/10.0, values[2], values[3:l]) # number, body, densityDesignVar, typeNumber, nodesindex):
                newElement.numberNodes = l - 3  # record the number of nodes
                #print g.nodes
                elementsList.append(newElement)

                count += 1
        return elementsList

    def writeElementMeshFile(self, rootFileDirectory,filename, elementsList):
        # ----------------------------
        # Export the new elements file.
        # ----------------------------
        f = os.path.join(rootFileDirectory,filename)
        with open(f, "wt") as fout:
            for element in elementsList:
                # print element.nodes
                nodeString = " ".join(map(str, element.nodesIndex))
                #print nodeString
                #print element.nodeNumber
                text = '{0} {1} {2} {3} \n'.format(element.nodeNumber, element.bodyNumber,504, nodeString)
                #print text
                fout.write(text)

    def readLocalStiffnessFile(self,path, filename, elementsList):
        f = os.path.join(path,filename)
        with open(f, "rt") as fin:
            count = 0
            countElement = -1  # the element indexes start at 0?? I think??

            currentElement = elementsList[0]
            newElement = re.compile('--')

            # loop over the lines of the file
            for line in fin:

                # if found a '--' then a new element
                m = newElement.search(line)
                if m:
                    found = m.group(0)
                    countElement = countElement + 1
                    # print countElement

                    # update the current element
                    currentElement = elementsList[countElement]
                    count = 0
                    continue  # go to the next line.

                # add the current value to the stiffness matrix
                currentElement.stiffnessMatrix.append(float(line))
                count = count + 1

    def readLocalMassFile(self,folderpath, filename, elementsList):
        f = os.path.join(folderpath,filename)
        with open(f, "rt") as fin:
            count = 0
            countElement = -1  # the element indexes start at 0?? I think??

            currentElement = elementsList[0]
            newElement = re.compile('--')

            # loop over the lines of the file
            for line in fin:

                # if found a '--' then a new element
                m = newElement.search(line)
                if m:
                    found = m.group(0)
                    countElement = countElement + 1
                    # print countElement

                    # update the current element
                    currentElement = elementsList[countElement]
                    count = 0
                    continue  # go to the next line.

                # add the current value to the mass matrix
                currentElement.massMatrix.append(float(line))
                count = count + 1


    def updateElementDensities(self,elementsList, nodeList, targetVolume):
        print 'ok'
        lower = 0.0
        upper = 1.0e1
        endCritiera = 1e-4
        move = 0.2
        curent_totalDensity = 0

        countloops = 0


        #while upper-lower > endCritiera:
        for i in range(1,20):
            middle = 0.5*(lower+upper)

            countloops+=1
            curent_totalDensity = 0.000
            totalMass = 0.00
            for e in elementsList:
                # designVar.w = max(0.001,max(designVar.w-move,min(1.,min(designVar.w+move,designVar.w.*sqrt(-g1./lmid)))));
                try:
                    t = e.densityDesignVar*math.sqrt(-e.sensitivity/middle)
                except:
                    t = 0
                    print 'sqrt error'
                print e.sensitivity,middle,math.sqrt(-e.sensitivity/middle),e.densityDesignVar,t
                newDensity = max(0.001,
                                  max(e.densityDesignVar-move,
                                      min(1,
                                          min(e.densityDesignVar + move,t

                                             )
                                          )
                                      )
                                )
                e.densityDesignVar = newDensity
                totalMass += newDensity
                # print e.densityDesignVar


            curent_totalDensity = totalMass/len(elementsList) # average the densities

            print 'loops ', countloops ,' middle ', middle, ' upper ', upper,' lower ', lower, ' cDensity ', curent_totalDensity,' targetD ', targetVolume
            if curent_totalDensity -targetVolume > 0:
                lower = middle
                print 1
            else:
                upper = middle

        # Update the body numbers for the elements.
        for e in elementsList:
            temp = e.densityDesignVar*9.0+2
            temp2 = min(int(round(temp)),10)
        e.bodyNumber = temp2 # 10 discritized bodies/materials
        #print temp, " ",  e.bodyNumber
        print 'some vars changed for debugging in updateElementDensities funciton, '
        print 'check the loop criteria, upper bounds, lower bounds, ....'


        #     l2 = lmid;
        #   end

        #
        # Matlab code of the optimal criteria
        #
        # %%%%%%%%%% OPTIMALITY CRITERIA UPDATE Gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # function [wnew]=OC_gradient(g1,designVar, settings)
        # l1 = 0; l2 = 100000; move = 0.05;
        # while (l2-l1 > 1e-4)
        #   lmid = 0.5*(l2+l1);
        #   designVar.w = max(0.001,max(designVar.w-move,min(1.,min(designVar.w+move,designVar.w.*sqrt(-g1./lmid)))));
        #
        # %   desvars = max(VOID, max((x - move), min(SOLID,  min((x + move),(x * (-dfc / lammid)**self.eta)**self.q))))
        #
        # [volume1, volume2] = designVar.CalculateVolumeFractions(settings);
        # currentTotal = volume1+volume2;
        # %currentvolume=volume1+volume2;
        #
        #   %if currentvolume - volfrac > 0;
        #   if currentTotal - settings.totalVolume > 0;
        #     l1 = lmid;
        #   else
        #     l2 = lmid;
        #   end
        #
        #   wnew =  designVar.w ;
        # end
