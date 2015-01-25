import math
import numpy as np

outerCircleRadius = 63.5 # mm
numsections = 11
degressPerSection = 360/float(numsections)

f = open('instance.xml','w')

for i in range(0,11):
    angle = degressPerSection*i+180
    angleRadians = math.radians(angle)
    x = math.cos(angleRadians)*outerCircleRadius
    y = math.sin(angleRadians)*outerCircleRadius
    print   'section: ' + str(i) + ' x = ' + str(x) + ' y= ' + str(y)
    
    f.write(' <instance objectid="{0}"> \n' .format(i))
    f.write('  <deltax>{0:.2f}</deltax> \n'.format(x))
    f.write('   <deltay>{0:.2f} </deltay>\n'.format(y))
    f.write('   <rz>0</rz>\n')
    f.write(' </instance>\n')
    
f.close()