step = 0.1
x = 0:step:5
y = 0:step:5
[x2,y2] = meshgrid(x,y)

z2 = sin(x2)+sin(y2)

surf2stl('test.stl',x2,y2,z2,'ascii');
