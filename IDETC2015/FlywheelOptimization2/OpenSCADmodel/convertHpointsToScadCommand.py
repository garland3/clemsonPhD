
import numpy as np
originalFile = 'hpoints1.csv'
outFile ='scadCommands2.txt'


with open(originalFile, "rt") as fin:
     # split the first line which is the x
     line = fin.readline().rstrip()
     x1 = line.split(',')
    
      # split the second line which is the y
     line = fin.readline().rstrip()
     y1 = line.split(',')


y = np.arange(0.02,0.2,0.001)
x = np.interp(y,x1,y1)   

     
with open(outFile, "wt") as fout:
     fout.write('polygon(points = [')
     count = 0 # keep track of how many points we add to the pointLocationList
     
     pointLocationList = ''
     pointsList = ''
     
     
     # start at the origin

     np.insert(x,0,'0.02')
     np.insert(x,0,'0')
     
     #x.insert(1,'0')
     #y.insert(1,'.02')
     
     # Go to the x axis at the end of the 
     np.append(x,'0.2')
     np.append(y,'0')
    
     
     # go to the bezier curve
     for pointX in x:
        pointLocationList+=('[' + str(pointX) + ',' + str(y[count]) + '],')   
        
        # Write the order that the points should be used. 
        pointsList+=(str(count) + ',')
        count=count+1      

      
    
     
     
     pointLocationList = pointLocationList[:-1] # strip the last comma
     fout.write(pointLocationList)
    
     
     pointsList = pointsList[:-1] # remove the    last comma
     fout.write('], paths=[[' + pointsList)
      
     fout.write(']]);')
    