import re
import numpy as np


class Config:
    def __init__(self):
        self.originalFile = 'block100InfillForGradient.gco'
        self.outFile = 'blockWithGradientStrip.gco'

        # X60.260 Y80.260
        self.OriginOfPart = ( 60.260, 80.260,0.260);

        # The first few lines of code are used to print around the perimeter to initialize the extrusion.
        self.linesOfInitializatin = 93

        self.CrossSectionArea = (1.75/2)*(1.75/2)*np.pi  # pi*r^2
        self.VolResponseDelay = 2.524468763  # based of the testing I did
        self.responseE =  self.VolResponseDelay / self.CrossSectionArea

        # mm step resolution of extruder motor steps
        self.stepE = .05



class GcodeLine:
    
    def __init__(self,text, count, config):
        self.Configuration = config
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
        # The total extrusion upto this section (sections are bracketed by G93 E0 commands)
        self.ESection  = previousGcodeLine.ESection

        # Esection plus the current E, by default use the previous gcodes Ecum
        self.ECumulative =previousGcodeLine.ECumulative
        self.ExtrudedVolume = 0
        self.X = previousGcodeLine.X
        self.Y = previousGcodeLine.Y
        self.Z = previousGcodeLine.Z

        # if it has a G1 command, then
        # http://stackoverflow.com/questions/4666973/how-to-extract-a-substring-from-inside-a-string-in-python
        # https://docs.python.org/2/howto/regex.html#regex-howto

        if self.text.find("G1") != -1:
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
        elif self.text.find('G92') != -1:
            if self.text.find("E") != -1:
                self.ESection= self.ECumulative


        xOld = previousGcodeLine.X
        yOld = previousGcodeLine.Y
        zOld = previousGcodeLine.Z
        xCur = self.X
        yCur = self.Y
        zCur = self.Z
        self.Distance = np.sqrt((xOld-xCur)*(xOld-xCur)+(yOld-yCur)*(yOld-yCur)+(zOld-zCur)*(zOld-zCur))

        if self.Distance!=0 and self.E !=0 :
                    oldE = previousGcodeLine.E
                    self.ExtrusionMultiplier = np.abs(oldE - self.E)/self.Distance

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
        self.ExtrudedVolume = self.Configuration.CrossSectionArea*self.ECumulative

def volumeFractionFunction(X,Y,Z):
    xDist = 159.740-60.260
    yDist = 129.740-80.2660

    X = X-60.260 # normalize to be at 0,0
    Y = Y - 80.260
    
    # if(np.mod(X,30) <=15):
    #    F = 100
    # else:
    #    F = 0
    F = (X+Y)*100/148.954

    
    return F


# The main executable for this program.
def main():
    count = 0;
    config = Config()

    GcodeList = []

    oldG = GcodeLine(';first gcode',0,config)
    oldG.SetValuesManually(config.OriginOfPart[0], config.OriginOfPart[1], config.OriginOfPart[2], 0, 0, 0, 0, 0, 0.05, 0)

    print '(g.ECumulative, g.ExtrudedVolume, g.ExtrusionMultiplier, g.text)'

    # Part 1, Recording the Gcodes
    with open(config.originalFile, "rt") as fin:
        for line in fin:
            count += 1
            if count <= config.linesOfInitializatin: # skip the first 107 lines
                continue

            g = GcodeLine(line, count,config)
            g.analyze(oldG)
            g.CalculateExtrudedVolumue()
            GcodeList.append(g) # store the gcode object in a list
            print (g.ECumulative, g.ExtrudedVolume, g.ExtrusionMultiplier, g.text)

            oldG = g

    maxE = GcodeList[-1].ECumulative
    # Part 2

    # make a list of E increments that
    Elist = np.arange(config.responseE,maxE, config.stepE )
    print 'Elist size'+str(Elist.size)



    XPointList = []
    YPointList = []
    ZPointList = []
    FractionPointList = []


    with open(config.outFile, "wt") as fout:
        # write the first 107 lines
        with open(config.originalFile, "rt") as fin:
            countInit = 1
            for line in fin:
                if countInit <= config.linesOfInitializatin:  # skip the first 107 lines
                    fout.write(line)
                countInit += 1


        count = 0
        count2 = 0
        for ETarget in Elist:
            code = '; extruded target = %.4f  \n' % ETarget
            fout.write(code)
            print ETarget
            for g1 in GcodeList[count:]:

                if g1.ECumulative >= ETarget:
                    d = g1.Distance

                    if count !=0:
                        Gprevious = GcodeList[count-1]
                    else:
                        Gprevious = GcodeList[count]

                    diffE = np.abs(ETarget-Gprevious.ECumulative)
                    distanceE = g1.ECumulative-Gprevious.ECumulative

                    #offSetE = diffVol/config.CrossSectionArea
                    #offSetDistance = offSetE/exMult

                    xList = [Gprevious.X, g1.X]
                    yList = [Gprevious.Y, g1.Y]

                    uList = [0, 1]
                    # interpolate between the two points using the parameter u
                    if d !=0:
                        newX = np.interp([diffE/distanceE], uList, xList)
                        newY = np.interp([diffE/distanceE], uList, yList)
                    else:
                        newX = Gprevious.X
                        newY = Gprevious.Y

                    volFraction = volumeFractionFunction(newX,newY,g1.Z)
                    XPointList.append(newX)
                    YPointList.append(newY)
                    ZPointList.append(g1.Z)
                    FractionPointList.append(volFraction)

                    code = 'G93 R%.1f \n' % volFraction
                    fout.write(code)

                    # since we wrote a G93 command, we need to write the next G1 command
                    for g2 in GcodeList[count2:]:
                        if g2.ECumulative >= (ETarget-config.responseE):
                            d = g2.Distance

                            # Get the prevous Gcode
                            if count2 !=0:
                                Gprevious = GcodeList[count2-1]
                            else:
                                Gprevious = GcodeList[count2]

                            diffE = np.abs((ETarget-config.responseE)-Gprevious.ECumulative)
                            distanceE = g2.ECumulative-Gprevious.ECumulative
                            #diffE = np.abs(ETarget-config.responseE)

                            xList = [Gprevious.X, g2.X]
                            yList = [Gprevious.Y, g2.Y]

                            uList = [0, 1]
                            # interpolate between the two points using the parameter u
                            if d !=0:
                                newX = np.interp([diffE/distanceE], uList, xList)
                                newY = np.interp([diffE/distanceE], uList, yList)
                            else:
                                newX = Gprevious.X
                                newY = Gprevious.Y



                            #Enew = (targetVol-config.VolResponseDelay)*exMult
                            if(d != 0):
                                code = 'G1 X%.4f Y%.4f Z%.4f E%.4f F%0.1f\n' % (newX,newY,g2.Z,(ETarget-config.responseE), g2.Feedrate)
                                fout.write(code)
                                break # break out of this loop becasue we are finished
                            else:
                                code = 'G1 E%.4f F%0.1f\n' % ((ETarget-config.responseE), g2.Feedrate)
                                fout.write(code)
                                break # break out of this loop becasue we are finished


                        else:
                             count2 +=1

                    break # break out of hte gcode searchingn look and go to the next volfact

                else:
                    count +=1


        # write the last 18 lines which are related to turnning everything off.
        for g in GcodeList[-18:]:
            fout.write(g.text)


    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from pylab import *


    fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #ax.scatter(XPointList, YPointList,FractionPointList)
    plt.scatter(XPointList, YPointList,c=FractionPointList, marker='o')
    plt.colorbar()
    
    

    plt.show()





if __name__ == "__main__":
    main()