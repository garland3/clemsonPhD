import os

path = os.path.dirname(os.path.realpath(__file__))
list = os.listdir(path)	

for item in list:
	print item