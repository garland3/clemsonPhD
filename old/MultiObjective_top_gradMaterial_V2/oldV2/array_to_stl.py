# from pylab import imread
from stl_tools import numpy2stl

# from scipy.misc import lena, imresize
# from scipy.ndimage import gaussian_filter
# import scipy

#A = 256 * imread("NASA.png")
# https://github.com/thearn/stl_tools

# A = 256 * imread("macro5.png")
# print A.size
# print A.shape

# A = A[:, :, 0]+ A[:, :, 1]+ A[:, :, 2] #+ 1.0*A[:,:, 0] # Compose RGBA channels to give depth
#A =   A[:, :, 1] #+ 1.0*A[:,:, 0] # Compose RGBA channels to give depth
#print A
# A = gaussian_filter(A, 1)  # smoothing
#A = A[:, :, 0]
# A = -A # inverse the colors
#A = A+10

# print 'max '+ str(scipy.amax(A))
# print  'min '+ str(scipy.amin(A))
# print  'average '+ str(scipy.average(A))
# print  'median '+ str(scipy.median(A))

# numpy2stl(A, "leafHole2.stl", scale=2, mask_val=80, solid=True, max_width=150, max_depth=150,max_height=2)

from numpy import genfromtxt
my_data = genfromtxt('completeStucture1.000000_macroIteration_1.csv', delimiter=',')
numpy2stl(my_data, "caketopper.stl", solid=True)
