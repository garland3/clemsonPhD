 % X and Y position of the control points for a 1 by 1 square

 clc
 close all
[xpoints,ypoints] = meshgrid(0:1/2:1,0:1/2:1);

% zpoints = xpoints.^2-ypoints.^2
%zpoints = [ 0	-0.01761	-1.5	-1.5	-0.12588	-1.0698	0	0.76488	1.48953	];
 zpoints = [ -0.43551	-0.3036	-0.05595	-0.3333	0.17751	0.01629	-0.8418	-0.68175	-6.17E-01];


zpoints_v2 = ones(3,3);
zpoints_v2(1,:) = zpoints(1:3);
zpoints_v2(2,:) = zpoints(4:6);
zpoints_v2(3,:) = zpoints(7:9);

doplot = 1;
[ydisplacment, cost,maxDisplacement]= FEALevelSet_2D_meso(xpoints,ypoints,zpoints_v2, doplot) % Call FEA for first time


