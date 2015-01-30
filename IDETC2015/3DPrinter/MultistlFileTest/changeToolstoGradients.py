import re
numlayers = 1

originalFile = 'halfbyhalfbyoneBox11.amf.gcode'
outFile ='halfbyhalfbyoneBoxV3_withGradient.gco'

originalFile = 'circular.amf.gcode'
outFile ='circular_gradient.gco'

# with open(outFile, "wt") as fout:
    # with open(originalFile, "rt") as fin:
        # for line in fin:
            # if (line.find("Z")  != -1):
                # numlayers+=1
                
# print numlayers

# ratio = float(100)/numlayers

# print ratio

# index = 1;

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