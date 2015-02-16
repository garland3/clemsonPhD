xpoints = 0:1:10; % X position of the control points


% Version 6 of the FEA, not incorporates a 2D level set function with
% control points. 
 ypoints = [0	0.82323	0.29682	0.03189	0.00636	-0.00435	-0.33795	-0.72885	-0.1689	-1.05174	-0.39495]

doplot = 1;
[ydisplacment, cost]= FEALevelSet(xpoints,ypoints,doplot); % Call FEA for first time

