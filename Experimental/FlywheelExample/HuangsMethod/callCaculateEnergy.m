function [objective] = callCaculateEnergy(inputvars)

innerR =0.02; % m,  inner shoudl be 0.02
outerR = 0.2; % m,  out should be 0.2
w = 630;%'omega'

numpoints = 5;

hpoints = inputvars(1,1:numpoints);
VolFractpoints =  inputvars(1,(numpoints+1):2*numpoints);

maxStress = CalculateEnergy(hpoints,VolFractpoints,w, innerR, outerR); % Joule, Newton*meters
objective = -maxStress; 