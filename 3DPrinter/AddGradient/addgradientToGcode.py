numlayers = 1

originalFile = '1.gco'
outFile ='earringsGradient.gco'

with open(outFile, "wt") as fout:
    with open(originalFile, "rt") as fin:
        for line in fin:
            if (line.find("Z")  != -1):
                numlayers+=1
                
print numlayers

ratio = float(100)/numlayers

print ratio

index = 1;

with open(outFile, "wt") as fout:
    with open(originalFile, "rt") as fin:
        for line in fin:
            if (line.find("Z")  != -1):
                index+=1
                fout.write(line)
                currentRatio = float(ratio*index)
                fout.write('G93 R%.1f \n' % currentRatio )
            else:
                fout.write(line)