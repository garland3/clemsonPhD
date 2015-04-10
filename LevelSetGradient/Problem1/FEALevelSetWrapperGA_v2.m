function [obj, cinequality] = FEALevelSetWrapperGA_v2(zpoints)

[xpoints,ypoints] = meshgrid(0:0.5:6,0:0.5:2); 
zpoints_v2 = xpoints;
[xdim, ydim] = size(zpoints_v2);


 %zpoints = [0.88564,0.87801,0.19888,0.99926,0.58114,0.9267,0.98697,0.98855,0.98715,0.99791,-0.00020904,-0.77478,-0.86254,-0.79403,-0.94306,-0.60567,-0.84454,0.95515,-0.039918,-0.39113,-0.77063,-0.54174,-0.42951,-0.64096,-0.50736,-0.528,-0.85935,0.49363,-0.46684,-0.95031,-0.83302,-0.55248,-0.45001,-0.17243,-0.44986,-0.44418,-0.075676,0.96908,-0.8129,-0.70698,-0.69448,-0.39489,-0.49969,-0.3296,0.99873,0.99998,0.99998,0.98501,-0.91432,-0.151,-0.79768,-0.23579,-0.62175,-0.18284,-0.45696];

zpoints = ones(1,xdim*ydim)*1.1;
doplot = 1;



count = 1;
for i = 1:xdim
    zpoints_v2(i,:) = zpoints(count:ydim+count-1);
    count = count + ydim;  

end




[ydisplacment, cost,maxDisplacement] = FEALevelSet_2D_v8(xpoints,ypoints,zpoints_v2, doplot); % Call FEA for first time

obj = cost
cinequality = maxDisplacement-0.2;  % max Displacemen is positive. The max displacement must be less than 0.1 inch
% cequality = ydisplacment+0.0623; % displacement will be negative, so add 0.0623 to make = 0
%feature('GetPid')





