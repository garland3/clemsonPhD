# Copyright Anthony Garland 2015
# Changes tool change Gcodes to gradients for the Big Builder 3D printer

import re
numlayers = 1

# Change these to references to your gcode and gcode output
originalFile = 'test.gcode'
outFile ='leaf3_withGradient.gco'



numExtruders = 10

with open(outFile, "wt") as fout:
    with open(originalFile, "rt") as fin:
        for line in fin:
            if (line.find("M104")  != -1): # prevent from removing the warm up commands by putting this at the beggining. 
                 fout.write(line)
            elif(line.find("M109")  != -1): # prevent from removing the warm up commands by putting this at the beggining. 
                 fout.write(line)
                 
            elif (line.find("T0")  != -1):
                currentRatio = 0
                fout.write('G93 R%.1f \n' % currentRatio )
            # elif(line.find("T1") != -1):
            elif(re.findall('\\bT1\\b', line) != []): # avoid comflict with finding T10 by using a regular expression search with boundaries
                currentRatio = float(100/numExtruders)*1
                fout.write('G93 R%.1f \n' % currentRatio )
            elif(line.find("T2") != -1):
                currentRatio = float(100/numExtruders)*2
                fout.write('G93 R%.1f \n' % currentRatio )
            elif(line.find("T3") != -1):
                currentRatio = float(100/numExtruders)*3
                fout.write('G93 R%.1f \n' % currentRatio )
            elif(line.find("T4") != -1):
                currentRatio = float(100/numExtruders)*4
                fout.write('G93 R%.1f \n' % currentRatio )
            elif(line.find("T5") != -1):
                currentRatio = float(100/numExtruders)*5
                fout.write('G93 R%.1f \n' % currentRatio )
            elif(line.find("T6") != -1):
                currentRatio = float(100/numExtruders)*6
                fout.write('G93 R%.1f \n' % currentRatio )                
            elif(line.find("T7") != -1):
                currentRatio = float(100/numExtruders)*7
                fout.write('G93 R%.1f \n' % currentRatio )
            elif(line.find("T8") != -1):
                currentRatio = float(100/numExtruders)*8
                fout.write('G93 R%.1f \n' % currentRatio )
            elif(line.find("T9") != -1):
                currentRatio = float(100/numExtruders)*9
                fout.write('G93 R%.1f \n' % currentRatio )
            elif(line.find("T10") != -1):
                currentRatio = float(100/numExtruders)*10
                print 'R%.1f \n' % currentRatio
                fout.write('G93 R%.1f \n' % currentRatio )
            else:
                fout.write(line)