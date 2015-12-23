#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# by Panos Mavrogiorgos, email: pmav99 <> gmail
 
from vtk import *
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
 
# The source file
file_name = "paraviewout0001.vtu"
 
# Read the source file.
reader = vtkXMLUnstructuredGridReader()
reader.SetFileName(file_name)
reader.Update() # Needed because of GetScalarRange


# Get the coordinates of nodes in the mesh
nodes_vtk_array= reader.GetOutput().GetPoints().GetData()

#Get the coordinates of the nodes and their temperatures
nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array)
x,y,z= nodes_nummpy_array[:,0] , nodes_nummpy_array[:,1] , nodes_nummpy_array[:,2]

#for i in x:
#	print i
#--------------
# Works!!!!!!
# ----------
print len(x)

	
#The "Displacement" field is the third scalar in my vtk file
displacement_vtk_array = reader.GetOutput().GetPointData().GetArray(0)

T = vtk_to_numpy(displacement_vtk_array)

print len(T)


count = 0
with open('nodesAndDisplacement.txt', "wt") as fout:
	for i in x:
		string = "{0} {1} {2} {3}\n".format(x[count], y[count],x[count],T[count])
		#print string
		count = count+1
		fout.write(string)
#for i in d_numpy_array:
#	print i
#--------------
# Works!!!!!!
# ----------
	
# #Draw contours
# npts = 100
# xmin, xmax = min(x), max(x)
# ymin, ymax = min(y), max(y)

# # define grid
# xi = np.linspace(xmin, xmax, npts)
# yi = np.linspace(ymin, ymax, npts)
# # grid the data
# Ti = griddata((x, y), T, (xi[None,:], yi[:,None]), method='cubic')  

# ## CONTOUR: draws the boundaries of the isosurfaces
# CS = plt.contour(xi,yi,Ti,10,linewidths=3) 

# ## CONTOUR ANNOTATION: puts a value label
# plt.clabel(CS, inline=1,inline_spacing= 3, fontsize=12, colors='k', use_clabeltext=1)

# plt.colorbar() 
# plt.show() 