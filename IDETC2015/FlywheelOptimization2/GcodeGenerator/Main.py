__author__ = 'Anthony G'
import math
import numpy as np
from scipy.interpolate import interp1d

# generate a generic loop function
# inputs are radius from center, volfraction composition, height
# call it Loop(x,x,x)


# Find the tallest point on the flywheel
# find the inner R and outer R
# loop from layer 1 to max, where max reaches the tallest incrememnt
# --- loop from inner to outer radius
# --- Test if the height at each radius is below the given set of input points
#-------- YES, then 
#----------------- Calculate the correct volume fraction based on the radius
# ----------------  call loop(x,x,x)
#---------NO, then keep going and do nothing.

#extrusionMultiplier = 0.053425 # 0.045

# 0.053425 was too much
# 0.045 worked, but still too much
# 0.39 was not enough
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
endRadius = 200*scale # 0.2 meters
extrudeThickness = .4479 # mm # 0.4479999999999933
extrudeThickness = .4 # mm # 0.4479999999999933

# ------------
# Z Parameters
# -----------
layerThickness = 0.20 # 0.2 mm
startZ = layerThickness*0.5 #
endZ = 100*scale # 0.01 meter

def depositeloop(radius,volfrac,height,fileHandle,radiuscount,layerCount):

    # rotate the start position with each radius. Every 5 mm radius, rotate pi
    rotationRaio = math.pi/5
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
    e = 20
    cmd+= 'T1\n'
    cmd+= 'G92 E0\n'
    cmd +=  'G1 X{0:.2f} Y{1:.2f} Z{2:.2f} E{3:.2f} F{4:.2f}\n'.format(10, 10,0, 0, ff)
    cmd +=  'G1 X{0:.2f} Y{1:.2f} Z{2:.2f} E{3:.2f} F{4:.2f}\n'.format(190, 10,0, e, ff)
    cmd +=  'G1 X{0:.2f} Y{1:.2f} Z{2:.2f} E{3:.2f} F{4:.2f}\n'.format(190, 180,0, e*2, ff)
    cmd +=  'G1 X{0:.2f} Y{1:.2f} Z{2:.2f} E{3:.2f} F{4:.2f}\n'.format(10, 180,0, e*3, ff)
    cmd +=  'G1 X{0:.2f} Y{1:.2f} Z{2:.2f} E{3:.2f} F{4:.2f}\n'.format(10, 10,0, e*4, ff)

    e = 15
    cmd+= 'T0\n'
    cmd+= 'G92 E0\n'
    cmd +=  'G1 X{0:.2f} Y{1:.2f} Z{2:.2f} E{3:.2f} F{4:.2f}\n'.format(11, 11,0, 0, ff)
    cmd +=  'G1 X{0:.2f} Y{1:.2f} Z{2:.2f} E{3:.2f} F{4:.2f}\n'.format(189, 11,0, e, ff)
    cmd +=  'G1 X{0:.2f} Y{1:.2f} Z{2:.2f} E{3:.2f} F{4:.2f}\n'.format(189, 179,0, e*2, ff)
    cmd +=  'G1 X{0:.2f} Y{1:.2f} Z{2:.2f} E{3:.2f} F{4:.2f}\n'.format(11, 179,0, e*3, ff)
    cmd +=  'G1 X{0:.2f} Y{1:.2f} Z{2:.2f} E{3:.2f} F{4:.2f}\n'.format(11, 11,0, e*4, ff)
    cmd+= 'T0\n'


    cmd+=depositeloop(endRadius +4,100,startZ,f,0,0)
    cmd+=depositeloop(endRadius +3.75,100,startZ,f,1,0)
    cmd+=depositeloop(endRadius +3.5,0,startZ,f,2,0)
    cmd+=depositeloop(endRadius +3.25,0,startZ,f,3,0)
    cmd+=depositeloop(endRadius +3,0,startZ,f,4,0)

    f.write(cmd)

def retract(amoutExtrude, filehandle):
    cmd = 'G1 E{0:.2f}\n'.format(amoutExtrude)
    filehandle.write(cmd)

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
    #cmd = moveToCenter()
    filehandle.write(intro)
   # filehandle.write(cmd)


# Read the z height points
originalFile = 'hpoints3.csv'
volFracOriginalFile = 'volFracpoints3.csv'

with open(originalFile, "rt") as fin:
    # split the first line which is the x (r)
    line = fin.readline().rstrip()
    r = [float(x) for x in line.split(',')]

    # split the second line which is the y (Z)
    line = fin.readline().rstrip()
    z = [float(y) for y in line.split(',')]

# read the vol fraction points
with open(volFracOriginalFile, "rt") as fin:
    # split the first line which is the x (r)
    line = fin.readline().rstrip()
    rVolFraction = [float(x) for x in line.split(',')]

    # split the second line which is the y (Z)
    line = fin.readline().rstrip()
    VolFraction = [float(y) for y in line.split(',')]


r = [r*1000 for r in r] # convert to mm from m
z = [z*1000 for z in z] # convert to mm from m

r = [r*scale for r in r] # convert to mm from m
z = [z*scale for z in z] # convert to mm from m


rVolFraction = [rVolFraction*1000 for rVolFraction in rVolFraction] # convert to mm from m
VolFraction = [VolFraction*100 for VolFraction in VolFraction] # convert from decimal to %

rVolFraction = [rVolFraction*scale for rVolFraction in rVolFraction] # convert to mm from m
#VolFraction = [VolFraction*scale for VolFraction in VolFraction] # convert to mm from m

print 'r'
print r
print 'z'
print z
print 'rVolFraction'
print rVolFraction
print 'VolFraction'
print VolFraction

# ---------------------------------
# Begin actually doing stuff
# --------------------------------






import matplotlib.pyplot as plt
#xnew = np.arange(.02, .2, 0.001)
#ynew = np.interp(xnew,x,y) # interp1d(x, y)
#f2 = interp1d(x, y, kind='cubic')

#print xnew
#plt.plot(x,y,'o',xnew,ynew,'-')  #, xnew)  , f2(xnew),'--')
#plt.legend(['data', 'linear'], loc='best')
#plt.legend(['data', 'linear', 'cubic'], loc='best')
#plt.show()

layerHeightList=np.arange(startZ,endZ,layerThickness)
radiusList = np.arange(startRadius,endRadius,extrudeThickness)
maxZatEachRadiusList = np.interp(radiusList,r,z)
volFractAtEachRadiusList = np.interp(radiusList,rVolFraction,VolFraction)

print 'layerHeight List {0}'.format( layerHeightList)
print 'Radius list {0}'.format(radiusList)
print 'Max Z at each Radius {0}'.format(maxZatEachRadiusList)
print 'VolFract at each radius {0}'.format(volFractAtEachRadiusList)

layerCount = 0
print 'Number of layers is ' + str(layerHeightList.size)

f = open('cylinder.gco','w')
writeBeginningCode(f)
cleanNozzle(f)
# retract(-1,f)

for zheight in layerHeightList:
    count = 0
    print layerCount
    changeZHeight(f,zheight)
    # retract(1,f)

    for i in radiusList:

        radius = i
        maxZallowed = maxZatEachRadiusList[count]
        localVolFract = volFractAtEachRadiusList[count]

        #volFract = i*10
        z = zheight
        #print 'Layer {0} , radius {1}, Max Allowed Z At Radius {2}, Current Z {3}\n'.format(layerCount,radius,maxZallowed,z)
        if(z<=maxZallowed):
            commands = depositeloop(radius,localVolFract,z,f,count,layerCount)
            f.write(commands)
        count+=1

    # retract(-1,f)
    layerCount+=1

endingCode(f)
f.close()

