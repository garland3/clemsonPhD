import re
import numpy as np

count = 0;
divisions = 10

#Material 1
#  Name = "Aluminium (generic)"
#  Heat Conductivity = 237.0
#  Youngs modulus = 70.0e9
#  Mesh Poisson ratio = 0.35
#  Heat Capacity = 897.0
#  Density = 2700.0
#  Poisson ratio = 0.35
#  Sound speed = 5000.0
#  Heat expansion Coefficient = 23.1e-6
#End

matValues = """Material {0}
  Name = "Aluminium (generic) {0}"
  Heat Conductivity = {1}
  Youngs modulus = {2}
  Mesh Poisson ratio = {3}
  Heat Capacity = {4}
  Density = {5}
  Poisson ratio = {6}
  Sound speed = {7}
  Heat expansion Coefficient = {8} \n"""
  
# Body 1
#  Target Bodies(1) = 1
#  Name = "Body 1"
#  Equation = 1
#  Material = 1
# End
  
bodyValues = """Body {0}
  Target Bodies(1) = {0}
  Name = "Body {0}"
  Equation = 1
  Material = {0}
End \n\n"""



max_Name = "Aluminium (generic)"
max_HeatConductivity = 237.0
max_Youngs_modulus = 70.0e9
max_MeshPoissonratio = 0.35
max_HeatCapacity = 897.0
max_Density = 2700.0
max_PoissonRatio = 0.35
max_SoundSpeed = 5000.0
max_HeatExpansionCoefficient = 23.1e-6

with open('material.txt', "wt") as fout:
	for i in range(0,divisions):
		#fout.write('Material ' + str(i+1) +'\n')
		
		fractionOfMax = i/float(divisions)
		#print fractionOfMax
		
		v=matValues.format(
			i+1, 
			fractionOfMax*max_HeatConductivity+1 ,
			fractionOfMax*max_Youngs_modulus+1 ,
			fractionOfMax*max_MeshPoissonratio+0.001,
			fractionOfMax*max_HeatCapacity+0.01,
			fractionOfMax*max_Density+0.01,
			fractionOfMax*max_PoissonRatio+0.0001,
			fractionOfMax*max_SoundSpeed+0.0001,
			fractionOfMax*max_HeatExpansionCoefficient+1e-9)
		fout.write(v)
		print v
		
		fout.write('End \n\n' )
		
		vv = bodyValues.format(
			i+1,
			i+1,
			i+1,
			i+1)
		print vv
		fout.write(vv)
		