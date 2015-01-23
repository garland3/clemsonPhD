
startRed=255
startGreen=192
startBlue=0

endRed=91
endGreen=155
endBlue=213


diffRed=endRed-startRed
diffGreen=endGreen-startGreen
diffBlue=endBlue-startBlue
print 'diff red: ' +  str(diffRed)
print 'diff green: ' +  str(diffGreen)
print 'diff blue: ' +  str(diffBlue)

numSections =11
div = numSections-1

for i in range(0,numSections):
    red = diffRed/float(div)*i +startRed
    green = diffGreen/float(div)*i + startGreen
    blue = diffBlue/float(div)*i + startBlue
    print 'section : '+str(i)+' red: '+str(red) + ' green: '+str(green)+ ' blue: ' +str(blue)