__author__ = 'Anthony G'



import os,sys
import Image
jpgfile = Image.open("picture.jpg")

print jpgfile.bits, jpgfile.size, jpgfile.format
