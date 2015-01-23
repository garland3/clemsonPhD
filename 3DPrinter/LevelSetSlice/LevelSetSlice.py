__author__ = 'Anthony G'
import numpy as np
import matplotlib.pyplot as plt

x = np.arange(-2, 2, 0.1)
y = np.arange(-2, 2, 0.1)
xx, yy = np.meshgrid(x, y, sparse=True)
# z = np.sin(xx**2 + yy**2) / (xx**2 + yy**2)
z = np.sin(xx)+np.sin(2*yy)

# normalize
maxZ =np.amax(z)
print maxZ
minZ = np.amin(z)
d = maxZ - minZ
z = z-minZ
Z = z/d


h = plt.contour(x,y,z)

print h.collections.count
print h.collections[0].get_paths()

plt._show()