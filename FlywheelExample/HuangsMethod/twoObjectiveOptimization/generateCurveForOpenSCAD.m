function generateCurveForOpenSCAD(hpoints,  iterationNum)
% plot(h(1,:),h(2,:));
filename = strcat('./GcodeGenerator/','hpoints',int2str(iterationNum),'.csv')
csvwrite(filename,hpoints);