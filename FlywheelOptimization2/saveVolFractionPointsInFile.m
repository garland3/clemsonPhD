function saveVolFractionPointsInFile(hpoints,  iterationNum)
% plot(h(1,:),h(2,:));
filename = strcat('./GcodeGenerator/','volFracpoints',int2str(iterationNum),'.csv')
csvwrite(filename,hpoints);