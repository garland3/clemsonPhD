function p = calculateDensity(volFracSn)
% table 3.2 in Huang's disertation, page 46
densitySN = 7.29e3; % kg/m^3
densityAL = 2.78e3; % kg/m^3
p = volFracSn*densitySN+(1-volFracSn)*densityAL;