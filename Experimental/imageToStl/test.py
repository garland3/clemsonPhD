#from stl_tools import numpy2stl

#from scipy.misc import lena, imresize
#from scipy.ndimage import gaussian_filter

import stl_tools 

import scipy.misc 
import scipy.ndimage 


A = scipy.imresize(lena(), (256, 256))  # load Lena image, shrink in half
A = gaussian_filter(A, 1)  # smoothing

numpy2stl(A, "examples/Lena.stl", scale=0.1, solid=False)