% Copyright Anthony Garland 2015
clear
clc
% xpoints = 0:1:10; % X position of the control points
[xpoints,ypoints] = meshgrid(0:2:10,0:0.5:1);

zpoints_v2 = ones(3,6);

zpoints = [ 0	-0.01761	-1.5	-1.5	-0.12588	-1.0698	0	0.76488	1.48953	-1.36128	-1.45662	-1.5	0	-0.34083	-1.49922	0.89559	-1.4997	-1.2336];

zpoints_v2(1,:) = zpoints(1:6);
zpoints_v2(2,:) = zpoints(7:12);
zpoints_v2(3,:) = zpoints(13:18);
doplot = 1;

% -------------------------------------------------
% - Make the voxels
% - Evaluted the phi function at each voxel (remember that the phi is only a
% function of xy, so apply it to all z voxels with the same xy)
% - Run the isosurface to get the vx
% - Run the stlwrite function to turn the vx into stl

% -----
% Make the voxels
isoCurves = -0.1:0.1:1
a = 4; % in
d = 7.5; % in
L = 10; % in
t = 0.5; % in 
h = 1; % in
step = 0.05;

% For testing make these really small
%L = 2; % in
%t = 0.5; % in 
%h = 1; % in
%step = 0.5;

xVox = 0:step:L;
yVox = 0:step:h;
zVox = 0:step:t;

[xVox2,yVox2,zVox2]=meshgrid(xVox,yVox,zVox);

[xVox2D,yVox2D] = meshgrid(xVox,yVox);

%tmpvol = zeros(20,20,20);       % Empty voxel volume
% splineZZ_v2 =  interp2(xpoints,ypoints,zpoints,splineXX_v2,splineYY_v2,'cubic'); % spline value at each x column
phiVox_zVox2D = interp2(xpoints,ypoints,zpoints_v2,xVox2D,yVox2D, 'linear');


bins = 10.0;

phiVox_zVox2D = round(phiVox_zVox2D*bins)/bins;




phiVox_zVox3D = xVox2; % get the same size
loopNum = size(phiVox_zVox3D,3); % Get the length of the Z dimmension

for i = 1:loopNum
    phiVox_zVox3D(:,:,i)=phiVox_zVox2D;
end

% phiVox_zVox3D(:,:,) = phiVox_zVox2D



count = 0;
for level = isoCurves
    fv = isosurface(xVox2,yVox2,zVox2,phiVox_zVox3D,level);  % Create the patch object
    isosurface(xVox2,yVox2,zVox2,phiVox_zVox3D,level);  % Create the patch object
    
    filename= strcat('results',num2str(count), '.stl'); % generate the filename
    %filename = char(filename)
   
    stlwrite(filename,fv, 'mode','ascii')        % Save to binary .stl
    
    count = count+1;
end