%hpoints = [0 1 2 3 4;0 0 0 0 0] % test data
%  hpoints = [0 1 2 3 4;0 0 0 0 0] % test data

csvwrite('hpoints.csv',hpoints)
csvwrite('VolFractpoints.csv',VolFractpoints)

innerR =0.02; % m,  inner shoudl be 0.02
outerR = 0.2; % m,  out should be 0.2

hpointsR = [innerR,0.04,0.1,0.12,outerR]; 
h = ones(2,5);
h(1,:) = hpointsR(1,:);
h(2,:) = hpoints(1,:);
h = transpose(h);

vf = ones(2,5);
vf(1,:) = hpointsR(1,:);
vf(2,:) = VolFractpoints(1,:);
vf = transpose(vf);


csvwrite('h.csv',h)
csvwrite('vf.csv',vf)


[mass, energy, volume, maxStress] = flywheel(h, vf)