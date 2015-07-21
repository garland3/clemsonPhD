__author__ = 'Anthony G'
# Copyright Anthony Garland 2015
# Modified to accept 10 control points

import math
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

endingCode = """G92 E0
M107
M104 S0
G28 X0
M84
M104 T0 S0
M140 S0
M84"""

# Adds a specific material volume fraction distribution in the z direction based off
# spline control points which are found through optimization.
# The spline is a function
# when it evaluates to 1, then it means 100% nylon
# when it evaluates to 0, then it means 100% PLA
# anything between is a mixture
# anything below 0 is empty

# The code will only modify an existing gcode file that prints a homogeneous beam
# the first 38 layers are the base ((PLA)
# layer 39 to 1308 is the beam that we will modify

volFracFile = 'hpoints_v5_beamdesign_11pnts.csv'
volFracFile = 'hpoints_v4_beamdesign.csv'
with open(volFracFile, "rt") as fin:
    # split the first line
    line = fin.readline().rstrip()
    volFractPointsY = [float(x) for x in line.split(',')]

print 'Volume Fraction Points Y'
print volFractPointsY

volFractPointsX =[0,2,4,6,8,10]
#volFractPointsX =[0,1,2,3,4,5,6,7,8,9,10]

xnew = np.arange(0, 10, 0.01)
ynew = np.interp(xnew,volFractPointsX,volFractPointsY) # interp1d(x, y)
f2 = interp1d(volFractPointsX, volFractPointsY, kind='cubic')

#print xnew
y2 = f2(xnew)

doplot = 1
if(doplot ==1):
    xLevelSet = [0, 5,10]
    yLevelSet = [0,0,0]

    #plt.plot(volFractPointsX,volFractPointsY,'o',xnew,ynew,'-', xnew,y2,'--',xLevelSet,yLevelSet,'-o')
    plt.figure( dpi = 600,figsize=(7,4))
   
    plt.plot(volFractPointsX,volFractPointsY,'o', xnew,y2,'--',xLevelSet,yLevelSet,'-o')
    plt.savefig('figOut.png', dpi = 600)
  
    # plt.legend(['data', 'linear'], loc='best')
    # plt.legend(['data', 'linear', 'cubic'], loc='best')
    #plt.legend(['data', 'cubic spline'], loc='best')
    # plt.show()


originalFile = 'v4_completeBeam.gco'
outFile ='v4_completeBeam_fgm.gco'


numBeamLayers = 1308-39 #  = 1269
beamLenth = 10
beamPerLayer =beamLenth/float(numBeamLayers)
print 'beamPerLayer'
print beamPerLayer

layerCount = 0
currentRatio = 100 # ratio 100 is 100% PLA on the back extruder
printLayer = 1 # set to one to print the layer, set to 0 to skip the layer

with open(outFile, "wt") as fout:
    with open(originalFile, "rt") as fin:
        for line in fin:
            if (line.find("Z")  != -1):

                # if we are in the first 38 layer, then we are printing the base. Do not change the currentRatio
                if(layerCount<39):
                    layerCount+=1
                    fout.write(line)
                    fout.write('G93 R%.1f \n' % currentRatio )
                else:
                    fout.write(line)
                    beamLayerCount = layerCount-39 # how many layers are on the actual beam
                    beamDistance = beamLayerCount*beamPerLayer # in the x direction
                    splineEvalutedAtDistance = f2(beamDistance)
                    #print 'beamLayerCount and splineEvalutedAtDistance'
                    print [beamLayerCount,beamDistance,np.array_str(splineEvalutedAtDistance) ]

                    # if below the level set
                    if(splineEvalutedAtDistance<0):
                        printLayer=0 # which will cause the layer to not print (we are inserting gcode)
                    else:
                        # if above the level set then specify the currect ratio
                        currentRatio = (1-splineEvalutedAtDistance)*100
                        fout.write('G93 R%.1f \n' % currentRatio )


                    layerCount+=1



            elif(line.find("T0") != -1 and printLayer==1):
                # Write the ratio code after a tool change, but do not increase the index because we are not changing z height.
                fout.write(line)
                fout.write('G93 R%.1f \n' % currentRatio )
            elif(printLayer==1):
                fout.write(line)
            # else:
                # skipped this gcode because we are not printing this layer

    # I removed the end of the code from the original gcode file and put it here to ensure it is
    # always in the new gcode
    fout.write(endingCode)
