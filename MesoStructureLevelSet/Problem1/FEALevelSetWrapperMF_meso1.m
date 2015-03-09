 % X and Y position of the control points for a 1 by 1 square
[xpoints,ypoints] = meshgrid(0:1/2:1,0:1/2:1);

% So we need to fill 16 values. a 4 by 4 square
zpoints_v2 = ones(3,3)
scale  = 1; 

zpoints = ones(1,9);


zpoints_v2(1,:) = zpoints(1:3)
zpoints_v2(2,:) = zpoints(4:6)
zpoints_v2(3,:) = zpoints(7:9)

doplot = 1;
[ydisplacment, cost,maxDisplacement]= FEALevelSet_2D_meso(xpoints,ypoints,zpoints_v2, doplot) % Call FEA for first time

