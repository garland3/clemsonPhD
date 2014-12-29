function generateCurveForOpenSCAD(hpoints,  iterationNum)
% plot(h(1,:),h(2,:));
filename = strcat('./OpenSCADmodel/','hpoints',int2str(iterationNum),'.csv')
csvwrite(filename,hpoints);