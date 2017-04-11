#!/usr/bin/evn python

import numpy as np
import scipy.linalg
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from numpy import genfromtxt    

# some 3-dim points
mean = np.array([0.0,0.0,0.0])
cov = np.array([[1.0,-0.5,0.8], [-0.5,1.1,0.0], [0.8,0.0,1.0]])
data = np.random.multivariate_normal(mean, cov, 50)

#xx=my_data = genfromtxt('../out1/ExxArrayForFitting1.csv', delimiter=',')
#yy=my_data = genfromtxt('../out1/EyyArrayForFitting1.csv', delimiter=',')
#zz=my_data = genfromtxt('../out1/RhoArrayForFitting1.csv', delimiter=',')


xx=my_data = genfromtxt('../out0/ExxArrayForFitting1.csv', delimiter=',')
#print xx
yy=my_data = genfromtxt('../out0/EyyArrayForFitting1.csv', delimiter=',')
#print yy
zz=my_data = genfromtxt('../out0/RhoArrayForFitting1.csv', delimiter=',')

data = np.c_[xx,yy,zz]
#print(data)


# https://gist.github.com/amroamroamro/1db8d69b4b65e8bc66a6


# regular grid covering the domain of the data
#X,Y = np.meshgrid(np.arange(-3.0, 3.0, 0.5), np.arange(-3.0, 3.0, 0.5))
mn = np.min(data, axis=0)
mx = np.max(data, axis=0)
X,Y = np.meshgrid(np.linspace(mn[0], mx[0], 20), np.linspace(mn[1], mx[1], 20))
XX = X.flatten()
YY = Y.flatten()

order = 2    # 1: linear, 2: quadratic, 3: cubit
if order == 1:
    # best-fit linear plane
    A = np.c_[data[:,0], data[:,1], np.ones(data.shape[0])]
    C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])    # coefficients
    
    # evaluate it on grid
    Z = C[0]*X + C[1]*Y + C[2]
    
    # or expressed using matrix/vector product
    #Z = np.dot(np.c_[XX, YY, np.ones(XX.shape)], C).reshape(X.shape)

elif order == 2:
    # best-fit quadratic curve
    A = np.c_[np.ones(data.shape[0]), data[:,:2], np.prod(data[:,:2], axis=1), data[:,:2]**2]
 
    C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])
    print C
    
    # evaluate it on a grid
    Z = np.dot(np.c_[np.ones(XX.shape), XX, YY, XX*YY, XX**2, YY**2], C).reshape(X.shape)
    
elif order == 3:
    # best-fit cubic curve
    A = np.c_[np.ones(data.shape[0]), data[:,:2], np.prod(data[:,:2], axis=1), data[:,:2]**2]
    print A
    C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])
    print C
    
    # evaluate it on a grid
    Z = np.dot(np.c_[np.ones(XX.shape), XX, YY, XX*YY, XX**2, YY**2], C).reshape(X.shape)
# plot points and fitted surface
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2)
ax.scatter(data[:,0], data[:,1], data[:,2], c='r', s=50)
plt.xlabel('X')
plt.ylabel('Y')
ax.set_zlabel('Z')
ax.axis('equal')
ax.axis('tight')
plt.show()