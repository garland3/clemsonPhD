

originalFile = 'hpoints1.csv'
outFile ='scadCommands.txt'


with open(originalFile, "rt") as fin:
     # split the first line which is the x
     line = fin.readline().rstrip()
     x = line.split(',')
    
      # split the second line which is the y
     line = fin.readline().rstrip()
     y = line.split(',')


with open(outFile, "wt") as fout:
     fout.write('polygon(points = [')
     count = 0 # keep track of how many points we add to the pointLocationList
     
     pointLocationList = ''
     pointsList = ''
     
     
     # start at the origin

     x.insert(0,'0.02')
     y.insert(0,'0')
     
     #x.insert(1,'0')
     #y.insert(1,'.02')
     
     # Go to the x axis at the end of the 
     x.append('0.2')
     y.append('0')
    
     
     # go to the bezier curve
     for pointX in x:
        pointLocationList+=('[' + pointX + ',' + y[count] + '],')   
        
        # Write the order that the points should be used. 
        pointsList+=(str(count) + ',')
        count=count+1      

      
    
     
     
     pointLocationList = pointLocationList[:-1] # strip the last comma
     fout.write(pointLocationList)
    
     
     pointsList = pointsList[:-1] # remove the    last comma
     fout.write('], paths=[[' + pointsList)
      
     fout.write(']]);')
    