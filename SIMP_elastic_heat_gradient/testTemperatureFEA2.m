% test temperature fea2

nelx = 10;
nely = 5;

x = ones(nely,nelx);
penal = 1;
temperatureFEA2(nelx,nely,x,penal)