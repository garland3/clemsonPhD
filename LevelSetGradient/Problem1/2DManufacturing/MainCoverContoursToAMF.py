__author__ = 'Anthony G'

# Copyright Anthony Garland 2015


import math
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from ContourPoint import ContourPoint

contourPointsFile = 'contours2D.csv'
with open(contourPointsFile, "rt") as fin:
    # split the first line
    line = fin.readline().rstrip()
    contourInfoX = [float(x) for x in line.split(',')]

    line = fin.readline().rstrip()
    contourInfoY = [float(x) for x in line.split(',')]


numPoints = len(contourInfoX)

count = 0
contourPointsRemaining = 0
isFirstPoint = False
zheight = -1

contourStartingPnts = []
contourEndingPnts = []

listOfContourPnts= []
contourNumber = 0

# http://www.mathworks.com/help/matlab/ref/contour-properties.html#prop_ContourMatrix
for x in contourInfoX:

    # if a new contour, then find out how many points it should have and the z height
    if(contourPointsRemaining ==0):
        contourNumber+=1
        print 'new contour'
        contourPointsRemaining = contourInfoY[count]
        zheight = x
        print 'zheight ' + str(zheight)
        isFirstPoint = True # add the next control point tot he list

    else:


        # if a contour point, then print it out.
        y = contourInfoY[count]
        print [x,y]

        # make a point and add it to the list
        pnt = ContourPoint(x,y,zheight,contourNumber)
        listOfContourPnts.append(pnt)

        if(isFirstPoint == True):
            # Then we are on the first point of this contour line
            contourStartingPnts.append(pnt)

            isFirstPoint = False

        # if on the last point of a contour line, then save the point
        if(contourPointsRemaining ==1):
            contourEndingPnts.append(pnt)

        contourPointsRemaining-=1
    count+=1
