from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

fig = plt.figure()
ax = fig.gca(projection='3d')
#X, Y, Z = axes3d.get_test_data(0.05)

size = 4
ZoffsetUp = 0

X = np.arange(-size, size, 0.25)
Y = np.arange(-size, size, 0.25)
X, Y = np.meshgrid(X, Y)
R = np.sqrt(X**2 + Y**2)
Z = np.sin(R)-0.4+X/4+Y/4+ZoffsetUp
Z2 = np.zeros_like(Z)
max = np.amax(Z)
min = np.amin(Z)
delta = max - min

Z = Z/delta*2
#Z = Z-(min -1)*2/delta
min2 = np.amin(Z)
print min2
Z = Z-min2 -1


V = 0

ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm)
cset = ax.contourf(X, Y, Z,V, zdir='z', offset=min, cmap=cm.coolwarm)

ax.contourf(X,Y,Z2)

ax.set_zlabel('Z')



ax.set_zlim(min, max)

plt.show()