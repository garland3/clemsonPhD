% Copyright Anthony Garland 2015 

targetDisplacment = -0.015;
maxNumIterations = 30;
moveLimit = 0.2;



% ----------------------------------------------------------------
% Start
% --------------------------------------------------------------

% original model
xpoints = 0:0.5:10; % X position of the control points
% ypoints = 1-xpoints.^2/80 % deceasing as you go right. Starts at 0.75

%ypoints = ones(1,21); %displacment is   -0.0104
 %ypoints = zeros(1,21);%    -0.0207
% ypoints = 0.5-0.5*xpoints.^2/25;
ypoints = 0+xpoints/10

[ydisplacment, strainEnergies,maxVonMisesInRegion,maxDisplacementInRegion]= FEALevelSet(xpoints,ypoints); % Call FEA for first time
figure(2) 
bar(strainEnergies)
strainEnergies

for i=1:maxNumIterations
    
    
end




