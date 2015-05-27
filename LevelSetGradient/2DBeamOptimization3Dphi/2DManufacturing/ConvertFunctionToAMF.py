__author__ = 'Anthony G'

import numpy as np
import matplotlib.pyplot as pl
import matplotlib.mlab as mlab

delta = 0.025
x = np.arange(-3.0, 3.0, delta)
y = np.arange(-2.0, 2.0, delta)
X, Y = np.meshgrid(x, y)
Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
# difference of Gaussians
Z = 10.0 * (Z2 - Z1)

cs = pl.contourf(X, Y, Z, 8)
# C = pl.contour(X, Y, f(X, Y), 8, colors='black', linewidth=.5)
#pl.show()


for collection in cs.collections:
    print collection
    for path in collection.get_paths():
        print path
        for polygon in path.to_polygons():
            print polygon.__class__
            print polygon
            pl.plot(np.array(polygon[:,0]),np.array(polygon[:,1]), 'x-')
            pl.show()

# "C:\Program Files (x86)\FreeCAD0.13\bin\FreeCADCmd.exe"
