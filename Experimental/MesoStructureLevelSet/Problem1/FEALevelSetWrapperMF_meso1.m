 % X and Y position of the control points for a 1 by 1 square
 % 
 % Version 1; 9 by 9 element unit cell with 9 control points 
 % Result:
 % zpoints = [ -0.43551	-0.3036	-0.05595	-0.3333	0.17751	0.01629	-0.8418	-0.68175	-6.17E-01];
 % Version 2, Updated the unit cell to have more elements. 

 clc
 close all
[xpoints,ypoints] = meshgrid(0:1/2:1,0:1/2:1);

% zpoints = xpoints.^2-ypoints.^2

% zpoints = [0.1976720849	-0.2016304502	-0.4018151903	0.652224579	0.1838078922	0.4381301333	-0.3062668771	0.9449075058	-1.221910625]


zpoints_v2 = ones(3,3);
zpoints_v2(1,:) = zpoints(1:3);
zpoints_v2(2,:) = zpoints(4:6);
zpoints_v2(3,:) = zpoints(7:9);

doplot = 1;
[ydisplacment, cost,maxDisplacement]= FEALevelSet_2D_mesoV2(xpoints,ypoints,zpoints_v2, doplot) % Call FEA for first time


