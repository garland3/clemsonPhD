import os 
import csv
from numpy import genfromtxt
import numpy
 
objectString = """<object id="%i" type="model">
   
      <color>
         <r>0.752941</r>
         <g>0.752941</g>
         <b>0.752941</b>
      </color>
      <mesh>
         <vertices>
            <vertex>
               <coordinates>
                  <x>0</x>
                  <y>0</y>
                  <z>0</z>
               </coordinates>
            </vertex>
            <vertex>
               <coordinates>
                  <x>1</x>
                  <y>0</y>
                  <z>0</z>
               </coordinates>
            </vertex>
            <vertex>
               <coordinates>
                  <x>0</x>
                  <y>0</y>
                  <z>10</z>
               </coordinates>
            </vertex>
            <vertex>
               <coordinates>
                  <x>1</x>
                  <y>0</y>
                  <z>10</z>
               </coordinates>
            </vertex>
            <vertex>
               <coordinates>
                  <x>0</x>
                  <y>1</y>
                  <z>0</z>
               </coordinates>
            </vertex>
            <vertex>
               <coordinates>
                  <x>1</x>
                  <y>1</y>
                  <z>0</z>
               </coordinates>
            </vertex>
            <vertex>
               <coordinates>
                  <x>0</x>
                  <y>1</y>
                  <z>10</z>
               </coordinates>
            </vertex>
            <vertex>
               <coordinates>
                  <x>1</x>
                  <y>1</y>
                  <z>10</z>
               </coordinates>
            </vertex>
         </vertices>
          <volume >
			<metadata type="slic3r.extruder">%i</metadata>
            <triangle>
               <v1>2</v1>
               <v2>6</v2>
               <v3>0</v3>
            </triangle>
            <triangle>
               <v1>6</v1>
               <v2>4</v2>
               <v3>0</v3>
            </triangle>
            <triangle>
               <v1>1</v1>
               <v2>2</v2>
               <v3>0</v3>
            </triangle>
            <triangle>
               <v1>4</v1>
               <v2>1</v2>
               <v3>0</v3>
            </triangle>
            <triangle>
               <v1>3</v1>
               <v2>2</v2>
               <v3>1</v3>
            </triangle>
            <triangle>
               <v1>5</v1>
               <v2>3</v2>
               <v3>1</v3>
            </triangle>
            <triangle>
               <v1>4</v1>
               <v2>5</v2>
               <v3>1</v3>
            </triangle>
            <triangle>
               <v1>3</v1>
               <v2>7</v2>
               <v3>2</v3>
            </triangle>
            <triangle>
               <v1>7</v1>
               <v2>6</v2>
               <v3>2</v3>
            </triangle>
            <triangle>
               <v1>5</v1>
               <v2>7</v2>
               <v3>3</v3>
            </triangle>
            <triangle>
               <v1>6</v1>
               <v2>7</v2>
               <v3>4</v3>
            </triangle>
            <triangle>
               <v1>7</v1>
               <v2>5</v2>
               <v3>4</v3>
            </triangle>
         </volume>
      </mesh>
	    </object>
	  """
	

startAMFfileString = """<?xml version="1.0" encoding="UTF-8"?>
	<amf unit="millimeter" version="1.1" xml:lang="en">\n"""

objectInstanceString = """<instance objectid="%i">
		
         <deltax>%i</deltax>
         <deltay>%i</deltay>
         <deltaz>0.00</deltaz>
         <rx>0</rx>
         <ry>0</ry>
         <rz>0</rz>
      </instance>\n"""

def writeConsellationObject(scaledValue,x,y):
	roundedValue = numpy.round(scaledValue,1)
	roundedValue*=10
	roundedValue+=1 # so that we start with object 1 and go to 11
	print roundedValue
	out = objectInstanceString % (roundedValue, x,y	)
	return out
	
	

def writeConstellation	():
	outString = ""
	outString+=  '<constellation id="100">\n'
	
	
	my_data = genfromtxt('gradAndStuct160_80by80.csv', delimiter=',')
	arraySize =  my_data.shape
	for x in range(1,arraySize[0]):
		for y in range(1,arraySize[1]):
			currentValue = my_data[x,y]
			# print currentValue
			
			# if = 1, then void, so if greater than 1, ...
			if(currentValue > 1):
			
				scaledCurrentValue = (currentValue-2)/2
				outString+= writeConsellationObject(scaledCurrentValue,x,y)
				
				
				
	#print arraySize[0]
	
	# if the value is 1 then void. 
	# if 2 = material 1
	# if 4 = material 2
	# between 2 and 4 is a gradient of material 1 and 2
	
	
	
	
	outString+=  ' </constellation>\n'
	return outString
	
	
	
endAMFFile = "</amf>"
outputFile = "test.amf"
ff = open(outputFile,'w')
ff.write(startAMFfileString)
ff.write(writeConstellation())

# Write the cube objects with different extruders. 
for i in range(1,12):
	objectXML = objectString % (i, i)
	ff.write(objectXML)
	
ff.write(endAMFFile)
ff.close()