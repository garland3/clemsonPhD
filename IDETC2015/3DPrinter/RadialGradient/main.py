__author__ = 'Anthony G'
import math
import numpy as np

filenameOut = 'ringGradient2.gco'

# extrusionMultiplier = 0.045
extrusionMultiplier = 0.05
firstLayerMultiplier=1.5
ff = 3000 # 3000
percentSlowerFirstLayer = 0.25
XOffset = 100
YOffset = 100
scale = .15 # scale 15%

# -----------
# Radius  Parameters
# ----------
# in mm
startRadius = 100*scale # 0.2 meters
endRadius = 250*scale # 0.4 meters
clearnNozzleRadius = endRadius+5
extrudeThickness = 0.35 # mm # 0.4479999999999933
numRings = (endRadius - startRadius)/extrudeThickness
print 'Number of rings is :' + str(numRings)

# ------------
# Z Parameters
# -----------
layerThickness = 0.20 # 0.2 mm
startZ = layerThickness #
endZ = 30*scale # 0.01 meter






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






    if (radiuscount == 0):
        d = 2 # extrude extra on the first inner radius
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

    code = 'G1 X{0:.2f} Y{1:.2f} Z{2:.2f} E{3:.2f} F{4:.2f}\n'.format(xcurrent+XOffset, ycurrent+YOffset,height, 0, ff)
    commands+=code

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

    endRadius = clearnNozzleRadius
    #startZ = ZHeightStep
    cmd+=depositeloop(endRadius +4,100,startZ,f,0,0)
    cmd+=depositeloop(endRadius +3.75,100,startZ,f,1,0)
    cmd+=depositeloop(endRadius +3.5,0,startZ,f,2,0)
    cmd+=depositeloop(endRadius +3.25,0,startZ,f,3,0)
    cmd+=depositeloop(endRadius +3,0,startZ,f,4,0)

    f.write(cmd)


def writeBeginningCode(filehandle):
    intro = """
G21
M104 T0 S215
M104 T1 S215
G28
M106
M109 T0 S240
M109 T1 S240
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

def Main():
    f = open(filenameOut,'w')
    writeBeginningCode(f)
    cleanNozzle(f)

    radiusList = np.arange(startRadius,endRadius,extrudeThickness)
    print 'end Z is :' + str(endZ)
    layerHeightList=np.arange(startZ,endZ,layerThickness)

    numRadii = radiusList.size
    ratio = float(100)/numRadii
    print 'length of radius list' + str(radiusList.size)
    print radiusList

    layercount  = 0

    for z in layerHeightList:
        radiusCount = 0
        for r in radiusList:
             localVolFract = radiusCount*ratio
             commands = depositeloop(r,localVolFract,z,f,radiusCount,layercount)
             f.write(commands)
             radiusCount += 1

        layercount +=1


    endingCode(f)




Main()