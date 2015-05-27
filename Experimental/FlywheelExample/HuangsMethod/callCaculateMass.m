function [maxStress] = callCaculateMass(inputvars)

numpoints = 5;

innerR =0.02; % m,  inner shoudl be 0.02
outerR = 0.2; % m,  out should be 0.2

hpoints = inputvars(1,1:numpoints);
VolFractpoints =  inputvars(1,(numpoints+1):2*numpoints);

maxStress =CalculateMass(hpoints,VolFractpoints, innerR, outerR); % kg