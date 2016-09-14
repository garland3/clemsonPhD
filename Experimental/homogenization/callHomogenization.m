lx = 1;
ly = 1;

E=1000 % elastic mod
v =0.3 %poisons ratio

% E2=1000 % elastic mod
% v2 =0.3 %poisons ratio


lambda = E*v/((1+v)*(1-2*v));
mu = E/(2*(1+v));

phi = 90;


x = randi([1 2],5)
%homogenize(1,1,[.01 2],[0.02 4],90,x)

homogenize(lx, ly, [lambda lambda/1000], [mu mu/1000], phi, x)