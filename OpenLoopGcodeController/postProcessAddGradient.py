import re
import numpy as np

CrossSectionArea = (1.75/2)*(1.75/2)*np.pi  # pi*r^2
VolResponseDelay = 2.524468763 # based of the testing I did

linearStepResolution = .2  #  mm step resolution when trying to calculate the volume fraction at each line segment
VolStepRes = linearStepResolution*CrossSectionArea # the corresponding volume amount





class GcodeLine:
    
    def __init__(self,text, count):
        self.text = text
        self.IsG1 = 0
        self.E = 0
        self.ESection = 0
        self.lineNumber = count

        self.Feedrate = 0
        self.ECumulative =0
        self.ExtrudedVolume = 0
        self.X = 0
        self.Y = 0
        self.Z = 0

        self.ExtrusionMultiplier = 0
        self.Distance = 0

    def analyze(self, previousGcodeLine):

        
        # set everything to start out, then populate with values as we go. 
        self.IsG1 = 0
        self.ExtrudeRatio = 0
        self.E = 0
        self.ESection  = previousGcodeLine.ESection # The total extrusion upto this section (sections are bracketed by G93 E0 commands)
        self.ECumulative =previousGcodeLine.ECumulative  # Esection plus the current E, by default use the previous gcodes Ecum
        self.ExtrudedVolume = 0
        self.X = previousGcodeLine.X
        self.Y = previousGcodeLine.Y
        self.Z = previousGcodeLine.Z


        
        
        # if it has a G1 command, then
        # http://stackoverflow.com/questions/4666973/how-to-extract-a-substring-from-inside-a-string-in-python
        # https://docs.python.org/2/howto/regex.html#regex-howto

        if (self.text.find("G1")  != -1):
            self.IsG1 = 0
            
            # see if it has a E
            m = re.search('E(.+?)\\s', self.text)
            if m:
                found = m.group(1)
                self.E = float(found)
                self.ECumulative = self.E+self.ESection # add the current E to the base value for the section

            # see if it has a F
            m = re.search('F(.+?)\\s', self.text)
            if m:
                found = m.group(1)
                self.Feedrate = float(found)
            else:
                self.Feedrate = previousGcodeLine.Feedrate





            # Search for an X value
            m = re.search('X(.+?)\\s', self.text)
            if m:
                found = m.group(1)
                self.X = float(found)

            # Search for an Y value
            m = re.search('Y(.+?)\\s', self.text)
            if m:
                found = m.group(1)
                self.Y = float(found)

             # Search for an Z value
            m = re.search('Z(.+?)\\s', self.text)
            if m:
                found = m.group(1)
                self.Z = float(found)


        # IF we come to a extruder position command, the update the cummulative total to 0
        elif (self.text.find("G92")  != -1):
            if(self.text.find("E") !=-1):
                 self.ESection= self.ECumulative


        xOld = previousGcodeLine.X
        yOld = previousGcodeLine.Y
        zOld = previousGcodeLine.Z
        xCur = self.X
        yCur = self.Y
        zCur = self.Z
        self.Distance = np.sqrt((xOld-xCur)*(xOld-xCur)+(yOld-yCur)*(yOld-yCur)+(zOld-zCur)*(zOld-zCur))

        if(self.Distance!=0 and self.E !=0):
                    oldE = previousGcodeLine.E
                    self.ExtrusionMultiplier= np.abs(oldE - self.E)/self.Distance



    def SetValuesManually(self,X,Y,Z,E,isG1,ESection,ECumulative,ExtrudedVolume,Multiplier, distance ):
        self.IsG1 = isG1
        self.E = E
        self.ESection = ESection
        self.ECumulative = ECumulative
        self.ExtrudedVolume = ExtrudedVolume
        self.X = X
        self.Y = Y
        self.Z = Z

        self.ExtrusionMultiplier = Multiplier
        self.Distance = distance

                
            
        
    def CalculateExtrudedVolumue(self):
        self.ExtrudedVolume = CrossSectionArea*self.ECumulative

def volumeFractionFunction(X,Y,Z):
    xDist = 159.740-60.260
    yDist = 129.740-80.2660

    X = X-60.260 # normalize to be at 0,0
    Y = Y - 80.260

    F = (X+Y)*100/148.954
    return F





originalFile = 'block100InfillForGradient.gco'
outFile = 'blockWithGradient2.gco'

# X60.260 Y80.260
OriginOfPart = ( 60.260, 80.260,0.260);
linesOfInitializatin = 93 # The first 107 lines of code are used to print around the perimeter to initialize the extrusion.

count = 0;

GcodeList = []

oldG = GcodeLine(';first gcode',0)
oldG.SetValuesManually(OriginOfPart[0], OriginOfPart[1], OriginOfPart[2], 0, 0, 0, 0, 0, 0.05, 0)

print '(g.ECumulative, g.ExtrudedVolume, g.ExtrusionMultiplier, g.text)'

with open(originalFile, "rt") as fin:
    for line in fin:
        count += 1
        if count <= linesOfInitializatin: # skip the first 107 lines
            continue

        g = GcodeLine(line, count)
        g.analyze(oldG)
        g.CalculateExtrudedVolumue()
        GcodeList.append(g) # store the gcode object in a list
        print (g.ECumulative, g.ExtrudedVolume, g.ExtrusionMultiplier, g.text)

        oldG = g

totalVolExtruded = GcodeList[-1].ExtrudedVolume
# Part 2

# make a list of extruded volumes were we want to change the composition
ExtrusionIncrements = np.arange(VolResponseDelay,totalVolExtruded,VolStepRes)
bookmark1 =0
bookmark2 = 0
previousG = GcodeList[0]
previousG2 =  GcodeList[0]
count = 0

XPointList = []
YPointList = []
ZPointList = []
FractionPointList = []


count = 1
with open(outFile, "wt") as fout:

    # write the first 107 lines
    with open(originalFile, "rt") as fin:
        for line in fin:
            if count <= linesOfInitializatin:  # skip the first 107 lines
                fout.write(line)
            count += 1

    for targetVol in ExtrusionIncrements:
        code = '; extruded target Vol = %.1f  \n' % targetVol
        fout.write(code)
        print targetVol
        for g in GcodeList[bookmark1:]:

            if g.ExtrudedVolume >= targetVol and g.Distance != 0:
                d = g.Distance
                exMult = g.ExtrusionMultiplier # E/distance = multiplierE

                oldExtruVol = previousG.ExtrudedVolume
                diffVol = np.abs(targetVol-oldExtruVol)

                offSetE = diffVol/CrossSectionArea
                offSetDistance = offSetE/exMult

                xList = [previousG.X, g.X]
                yList = [previousG.Y, g.Y]

                uList = [0, 1]
                # interpolate between the two points using the parameter u
                if(d!=0 and offSetDistance!=0):
                    newX = np.interp([offSetDistance/d], uList, xList)
                    newY = np.interp([offSetDistance/d], uList, yList)


                    volFraction = volumeFractionFunction(newX,newY,g.Z)
                    XPointList.append(newX)
                    YPointList.append(newY)
                    ZPointList.append(g.Z)
                    FractionPointList.append(volFraction)

                    code = 'G93 R%.1f \n' % volFraction
                    fout.write(code )
                    #print code

                    # since we wrote a G93 command, we need to write the next G1 command
                    for g2 in GcodeList[bookmark2:]:
                        if(targetVol == VolResponseDelay):
                            break

                        if(g2.ExtrudedVolume>=(targetVol-VolResponseDelay)):
                            d = g2.Distance
                            exMult = g2.ExtrusionMultiplier

                            oldExtruVol = previousG2.ExtrudedVolume
                            diffVol = np.abs( targetVol-oldExtruVol)

                            offSetE = diffVol/CrossSectionArea
                            offSetDistance = offSetE/exMult

                            xList = [previousG2.X, g2.X]
                            yList = [previousG2.Y, g2.Y]

                            uList = [0, 1]
                            # interpolate between the two points using the parameter u
                            if(d!=0 and offSetDistance!=0):
                                newX = np.interp([offSetDistance/d], uList, xList)
                                newY = np.interp([offSetDistance/d], uList, yList)
                            else:
                                newX = previousG2.X
                                newY = previousG2.Y



                            Enew = (targetVol-VolResponseDelay)*exMult
                            if(Enew != 0):
                                code = 'G1 X%.4f Y%.4f Z%.4f E%.4f F%0.1f\n' % (newX,newY,g2.Z,Enew, g2.Feedrate)
                                fout.write(code)
                                break # break out of this loop becasue we are finished
                            else:
                                 previousG2 = g2
                                 bookmark2+=1


                        else:
                             previousG2 = g2
                             bookmark2+=1

                    break # break out of hte gcode searchingn look and go to the next volfact


                else:
                    previousG = g
                    bookmark1+=1
            else:
                previousG = g
                bookmark1+=1

    # write the last 18 lines which are related to turnning everything off.
    for g in GcodeList[-18:]:
        fout.write(g.text)


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(XPointList, YPointList,FractionPointList)

#plt.show()



