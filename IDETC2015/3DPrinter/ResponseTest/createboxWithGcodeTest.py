# Copyright Anthony Garland 2015
import math
import numpy as np

filenameOut = 'responseTestPlate2.gco'

maxZ= 5 
theta = 45 # rotate the extruded rectangle this number of degrees. 
thetaR = math.radians(theta) # convert to radians
boxWidth  = 10
boxLength = 100
xStart = 100 -boxWidth/2
yStart = 100-boxLength/2
XOffset = 100
YOffset = 100

clearnNozzleRadius = math.sqrt(boxWidth*boxWidth*0.25+boxLength*boxLength*0.25)
print clearnNozzleRadius

extrusionMultiplier = 0.045 # 0.045
firstLayerMultiplier = 1.5
ff = 3000 # 3000
percentSlowerFirstLayer = 0.25

filamentWidth = 0.5
ZHeightStep = 0.24

numberofLayers = int(math.ceil(maxZ/ZHeightStep))
numberofLayers = 2


endZ = numberofLayers*ZHeightStep

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
    cmd = ''
    cmd+= 'T0\n'

    endRadius = clearnNozzleRadius
    startZ = ZHeightStep
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
    
     #Draw a whole zigzag
    z=0.24
    d = 0
    
    cmd= 'T0 \n G92 E0\n' # reset
    cmd += 'G1 X{0:.2f} Y{1:.2f} Z{2:.2f} E{3:.2f} F{4:.2f}\n'.format(xStart, yStart, z, 0, 5000)      # go to the start
    cmd += 'G1  E{0:.2f} F{1:.2f}\n'.format( 2, 1200)      # extrude 2 mm to start
    cmd+= 'G92 E0\n' # reset
    f.write(cmd)

    xOld = 78.898 
    yOld = 63.047

    filamentWidth = 2

    numLoop = int(boxWidth/(2*filamentWidth))
    actualWidth = (numLoop-1)*2*filamentWidth+filamentWidth


    #ratio = float(99)/numLoop
    
    codeVolFraction = 'G93 R{0:.2f} \n'.format(0)
    f.write(codeVolFraction)

    #print 'ratio per loop is ' + str(ratio)

    print(numLoop)
    print (actualWidth)
    print ("\n\n")

    for zz in range(1,numberofLayers):
        print z
        
        for j in range(0,numLoop):
            
            #   volfrac = float(ratio*j)
          
            if(j>=numLoop/2):            
                codeVolFraction = 'G93 R{0:.2f} \n'.format(100)
                f.write(codeVolFraction)
            
            if (zz ==1):
                extrusionMultiplierLocal = extrusionMultiplier*firstLayerMultiplier
                ffLocal = ff*percentSlowerFirstLayer
            else:
                extrusionMultiplierLocal = extrusionMultiplier
                ffLocal = ff
            
            print 'Extrusion multiplier ' + str(extrusionMultiplierLocal)
            x = 2*j*filamentWidth+xStart
            y = yStart
            d+= math.sqrt((x-xOld)*(x-xOld)+(y-yOld)*(y-yOld))*extrusionMultiplierLocal
            code = 'G1 X{0:.2f} Y{1:.2f} Z{2:.2f} E{3:.2f} F{4:.2f}\n'.format(x, y, z, d, ffLocal)          
            xOld = x
            yOld = y
            f.write(code)
            
            
            y = yStart+boxLength
            d+= math.sqrt((x-xOld)*(x-xOld)+(y-yOld)*(y-yOld))*extrusionMultiplierLocal
            code = 'G1 X{0:.2f} Y{1:.2f} Z{2:.2f} E{3:.2f} F{4:.2f}\n'.format(x, y, z, d,ffLocal)       
            xOld = x
            yOld = y
            f.write(code)
            
            x = (2*j+1)*filamentWidth+xStart          
            d+= math.sqrt((x-xOld)*(x-xOld)+(y-yOld)*(y-yOld))*extrusionMultiplierLocal
            code = 'G1 X{0:.2f} Y{1:.2f} Z{2:.2f} E{3:.2f} F{4:.2f}\n'.format(x, y, z, d,ffLocal)    
            xOld = x
            yOld = y
            f.write(code)
            
            
            y = yStart
            d+= math.sqrt((x-xOld)*(x-xOld)+(y-yOld)*(y-yOld))*extrusionMultiplierLocal
            code = 'G1 X{0:.2f} Y{1:.2f} Z{2:.2f} E{3:.2f} F{4:.2f}\n'.format(x, y, z, d,ffLocal)    
            xOld = x
            yOld = y
            f.write(code)
            f.write('\n')
        z+=ZHeightStep
        f.write('\n')
        

    # f.write('\n\n\n\n\n')
    # # Write the sides of the box
    # for zz in range(0,60):
        # z += ZHeightStep
        
       
        # x = xStart
        # y = yStart
        # d+= math.sqrt((x-xOld)*(x-xOld)+(y-yOld)*(y-yOld))*extrusionMultiplier
        # code = 'G1 X{0:.2f} Y{1:.2f} Z{2:.2f} E{3:.2f} F{4:.2f}\n'.format(x, y, z, d,ff)            
        # xOld = x
        # yOld = y
        # f.write(code)
        
        
        # y = yStart+boxLength
        # d+= math.sqrt((x-xOld)*(x-xOld)+(y-yOld)*(y-yOld))*extrusionMultiplier
        # code = 'G1 X{0:.2f} Y{1:.2f} Z{2:.2f} E{3:.2f} F{4:.2f}\n'.format(x, y, z, d,ff)    
        # xOld = x
        # yOld = y
        # f.write(code)
        
        # x = actualWidth +xStart
        # d+= math.sqrt((x-xOld)*(x-xOld)+(y-yOld)*(y-yOld))*extrusionMultiplier
        # code = 'G1 X{0:.2f} Y{1:.2f} Z{2:.2f} E{3:.2f} F{4:.2f}\n'.format(x, y, z, d,ff)    
        # xOld = x
        # yOld = y
        # f.write(code)
        
        
        # y = yStart
        # d+= math.sqrt((x-xOld)*(x-xOld)+(y-yOld)*(y-yOld))*extrusionMultiplier
        # code = 'G1 X{0:.2f} Y{1:.2f} Z{2:.2f} E{3:.2f} F{4:.2f}\n'.format(x, y, z, d,ff)     
        # xOld = x
        # yOld = y
        # f.write(code)

    endingCode(f)

Main()