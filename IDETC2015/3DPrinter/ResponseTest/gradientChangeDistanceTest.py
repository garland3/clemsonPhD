# Copyright Anthony Garland 2015

__author__ = 'Anthony G'
import math
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

zigzagDistance = 60

extrusionMultiplier = 0.045
firstLayerMultiplier=1.5
ff = 1200 # 3000
XOffset = 100
YOffset = 100
scale = .15 # scale 15%

# -----------
# Radius  Parameters
# ----------
# in mm
startRadius = 20*scale # 0.02 meters
endRadius = zigzagDistance*math.sqrt(2)+1
extrudeThickness = .4479 # mm # 0.4479999999999933
extrudeThickness = .4 # mm # 0.4479999999999933

# ------------
# Z Parameters
# -----------
layerThickness = 0.20 # 0.2 mm
startZ = layerThickness*0.7 #
endZ = 100*scale # 0.01 meter




def depositeLines(listOfXlocations, listofYLocations, listofZlocations, listOfVolFraction, fileHandle, listOfFF, layerCount):
    count = 0
    xOld = 0
    yOld = 0
    d= 0
    commands = ''
    
    
    extrusionMultiplierLocal = extrusionMultiplier
    if(layerCount == 0):
        extrusionMultiplierLocal = extrusionMultiplier*firstLayerMultiplier
    
    for xvalu in listOfXlocations:
        xcurrent = xvalu
        ycurrent =listofYLocations[count]
        height = listofZlocations[count]
        ff = listOfFF[count]
        volfrac=listOfVolFraction[count]
        
        codeVolFraction = 'G93 R{0:.2f} \n'.format(volfrac)
        commands+=codeVolFraction

        d+= math.sqrt((xcurrent-xOld)*(xcurrent-xOld)+(ycurrent-yOld)*(ycurrent-yOld))*extrusionMultiplierLocal        
        code = 'G1 X{0:.2f} Y{1:.2f} Z{2:.2f} E{3:.2f} F{4:.2f}\n'.format(xcurrent+XOffset, ycurrent+YOffset,height, d, ff)
        commands+=code
        xOld = xcurrent
        yOld = ycurrent


        count+=1
    return commands    


def depositeloop(radius,volfrac,height,fileHandle,radiuscount,layerCount):

    # rotate the start position with each radius. Every 5 mm radius, rotate pi
    rotationRaio = math.pi/5
    rotationRaio = 0 # set to zero so we can see a clean location of where 1 loop starts and the next one also starts for measuring. 
    tStart = 0+radius *rotationRaio
    tEnd = 2*math.pi+radius *rotationRaio

    stepsize = 0.05 # 0.2
    t = np.arange(tStart,tEnd,stepsize)
    t = np.append (t,tEnd) # make sure we get the last 2*pi in.
    x = np.cos(t)*radius
    y = np.sin(t)*radius

    commands = ''

    # Reset the extrusion distance
    codeResetExtrussion = 'G92 E0 \n'
    commands+=codeResetExtrussion

    # Write the vol fraction command
    codeVolFraction = 'G93 R{0:.2f} \n'.format(volfrac)
    commands+=codeVolFraction

    # Write the actual circle commands
    count = 0
    xOld = 0
    yOld = 0


    # move to the start of the circle with a g2 commammand
    xcurrent =np.cos(tStart)*radius
    ycurrent =np. sin(tStart)*radius


    ff = 5000
    code = 'G1 X{0:.2f} Y{1:.2f} Z{2:.2f} E{3:.2f} F{4:.2f}\n'.format(xcurrent+XOffset, ycurrent+YOffset,height, 0, ff)

    commands+=code
    if (radiuscount == 0):
        d = 4 # extrude extra on the first inner radius
    else:
        d = 0

    extrusionMultiplierLocal = extrusionMultiplier
    if(layerCount == 0):
        extrusionMultiplierLocal = extrusionMultiplier*firstLayerMultiplier


    xOld = xcurrent
    yOld = ycurrent

    if(radius <=4):
        ff = 1000
    elif(radius<8):
        ff = 2000
    else :
        ff = 6000

    for xvalu in x:
        xcurrent = xvalu
        ycurrent = y[count]

        d+= math.sqrt((xcurrent-xOld)*(xcurrent-xOld)+(ycurrent-yOld)*(ycurrent-yOld))*extrusionMultiplierLocal
        #commands+= 'G1 X' + str(xcurrent) + ' Y' +  str(ycurrent) + ' \n'
        code = 'G1 X{0:.2f} Y{1:.2f} Z{2:.2f} E{3:.2f} F{4:.2f}\n'.format(xcurrent+XOffset, ycurrent+YOffset,height, d, ff)
        commands+=code
        xOld = xcurrent
        yOld = ycurrent


        count+=1
    return commands

    print commands
    #print t
    #print x
    #print y

def cleanNozzle(f):
    ff =3000
    cmd = '\n'   


    cmd+=depositeloop(endRadius +4,100,startZ,f,0,0)
    cmd+=depositeloop(endRadius +3.75,100,startZ,f,1,0)
    cmd+=depositeloop(endRadius +3.5,0,startZ,f,2,0)
    cmd+=depositeloop(endRadius +3.25,0,startZ,f,3,0)
    cmd+=depositeloop(endRadius +3,0,startZ,f,4,0)

    f.write(cmd)



def changeZHeight(file, newZ):
    cmd = 'G1 Z{0:.2f}\n'.format(newZ)
    f.write(cmd)

def moveToCenter():
    cmd = 'G92 E0 \n'
    cmd+='G1 X100 Y100 F12000 \n'
    cmd+='G92 X0 Y0 E0 \n'
    return cmd

def endingCode(fileHandle):
    cmdEnd = """
;End GCode
M104 S0                     ;extruder heater off
M140 S0                     ;heated bed heater off (if you have it)
G91                                    ;relative positioning
G1 E-1 F300                            ;retract the filament a bit before lifting the nozzle, to release some of the pressure
G1 Z+0.5 E-5 X-20 Y-20 F9000 ;move Z up a bit and retract filament even more
G28 X0 Y0                              ;move X/Y to min endstops, so the head is out of the way
M84                         ;steppers off
G90                         ;absolute positioning """
    cmd  = "G0 F9000 X82.82 Y85.67 Z{0}\n".format(endZ+10)

    fileHandle.write(cmd)
    fileHandle.write(cmdEnd)


def writeBeginningCode(filehandle):
    intro = """
G21
M104 T0 S215
M104 T1 S215
G28
M106
M109 T0 S215
M109 T1 S215
G90
G92 E0
M82
T0
G92 E0
G1 E-4 F780
G92 E0
G1 Z0.26 F4200
 """

    filehandle.write(intro) 








layerCount = 0


f = open('responseTestCircle.gco','w')
writeBeginningCode(f)
cleanNozzle(f)

ff = 5000
    
xcurrent =0
ycurrent=0
d= 0
height=startZ
code = 'G1 X{0:.2f} Y{1:.2f} Z{2:.2f} E{3:.2f} F{4:.2f}\n'.format(xcurrent+XOffset, ycurrent+YOffset,height, d, ff)
volfrac=0
code+=  'G93 R{0:.2f} \n'.format(volfrac)
endRadius = endRadius/math.sqrt(2)

# listOfXlocations = [endRadius,-endRadius,endRadius,-endRadius,]
# listofYLocations = [endRadius,-endRadius,-endRadius,endRadius]
# listofZlocations = [startZ,startZ,startZ,startZ]
# listOfVolFraction = [0,100,100,0]
# listOfff = [1000,1000,1000,1000]

listOfXlocations = [0]
listofYLocations = [0]
listofZlocations = [startZ]
listOfVolFraction = [0]
listOfff = [2000]

distanceBetweenZags = 2
yD = int(math.ceil(zigzagDistance/distanceBetweenZags))

for t in range(-yD,yD):
    
    if(t % 2 ==0):
        listOfXlocations.append(zigzagDistance)
        listOfVolFraction.append(0)
    else:
         listOfXlocations.append(-zigzagDistance)
         listOfVolFraction.append(100)
         
    listofYLocations.append(t*distanceBetweenZags)
    listofZlocations.append(startZ)
    listOfff.append(2000)
    
    
    
    


code+= 'G92 E0 \n'

# code+=depositeLines(listOfXlocations, listofYLocations, listofZlocations, listOfVolFraction, f, listOfff, layerCount)

depositingRadius = int(endRadius-2)
print 'depositing radius is '+  str(depositingRadius)


code+=depositeloop(depositingRadius ,100,startZ,f,0,0)
code+=depositeloop(depositingRadius -1,100,startZ,f,1,0)
code+=depositeloop(depositingRadius -2,0,startZ,f,2,0)
code+=depositeloop(depositingRadius -3,0,startZ,f,3,0)
code+=depositeloop(depositingRadius -4,0,startZ,f,4,0)

f.write(code)

endingCode(f)
f.close()

