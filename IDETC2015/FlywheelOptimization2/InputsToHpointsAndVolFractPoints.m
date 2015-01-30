function [hpoints,VolFractpoints] = InputsToHpointsAndVolFractPoints(x)


innerR =0.02; % m,  inner shoudl be 0.02
outerR = 0.2; % m,  out should be 0.2


% hpoints will be the first 10
% volfractpoints will be the last 10
numberInputs = size(x,2);
numpointsForProfile = numberInputs/2;

% make the r distances in the radial direction. Evenly space them. 
length = outerR-innerR;
hpointsR = innerR:length/(numpointsForProfile-1):outerR;


hpoints = ones(2,numpointsForProfile);
hpoints(1,:) = hpointsR(1,:); % first row are the radial direction coordinate values
hpoints(2,:) = x(1,1:numpointsForProfile); % second row, are the height direction values
hpoints = transpose(hpoints); % transpose because all the equations expect a vertical list


VolFractpoints = ones(2,numpointsForProfile);
VolFractpoints(1,:) = hpointsR(1,:); % place the vol fract points at the same radial direction locations as the height functions
VolFractpoints(2,:) = x(1,(numpointsForProfile+1):2*numpointsForProfile);
VolFractpoints = transpose(VolFractpoints);

