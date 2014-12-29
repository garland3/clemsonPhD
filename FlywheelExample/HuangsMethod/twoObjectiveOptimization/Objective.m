function objective = Objective(input,w1,w2)

% assume 'input' has 20 values, 

innerR =0.02; % m,  inner shoudl be 0.02
outerR = 0.2; % m,  out should be 0.2
w = 630;%'omega'
shouldplot = 0; % set to 1 to see plots

[hpoints,VolFractpoints] = InputsToHpointsAndVolFractPoints(input);

s = -CalculateEnergy(hpoints, VolFractpoints,w,innerR, outerR);
p = CalculateMaxStress(hpoints, VolFractpoints, shouldplot);


objective = w1*(s+248277/198277)+w2*(p-18.7188e6)/(45e6-18.7188e6);


