# Copyright Anthony Garland 2015
# Adds a gradient in the Z direction existing code. 


# Change these two lines in order to add the gradient
originalFile = 'cylinder.gco'
outFile ='cylinder_gradient.gco'


numlayers = 1
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
            elif(line.find("T0") != -1):
                # Write the ratio code after a tool change, but do not increase the index because we are not changing z height. 
                fout.write(line)
                currentRatio = float(ratio*index)
                fout.write('G93 R%.1f \n' % currentRatio )
            else:
                fout.write(line)