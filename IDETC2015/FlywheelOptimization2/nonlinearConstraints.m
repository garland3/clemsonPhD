function [c,ceq] = nonlinearConstraints(x)
%c = ...     % Compute nonlinear inequalities at x.
% ceq = ...   % Compute nonlinear equalities at x.
ceq = [];
minAllowableEnergy = 50e3;
maxAllowableStress = 45e6; % 45MPa
maxAllowableMass = 75; 



innerR =0.02; % m,  inner shoudl be 0.02
outerR = 0.2; % m,  out should be 0.2
w = 630;%'omega'
shouldplot = 0; % set to 1 to see plots
[hpoints,VolFractpoints] = InputsToHpointsAndVolFractPoints(x);

c(1) = -(CalculateEnergy(hpoints, VolFractpoints,w,innerR, outerR)-minAllowableEnergy); % <=0
c(2) = CalculateMaxStress(hpoints, VolFractpoints, shouldplot)-maxAllowableStress;% <=0
c(3) = CalculateMass(hpoints, VolFractpoints,innerR, outerR) - maxAllowableMass; % <=0

