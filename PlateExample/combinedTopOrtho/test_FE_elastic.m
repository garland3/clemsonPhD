nelx = 20;
nely = 10;
x = ones(nely,nelx);
F = sparse(2*(nely+1)*(nelx+1),1);
penal=1;


nodeNumber = 2*(nely+1)*floor( (nelx+1)/2);
F(nodeNumber,1) = -10; % force at particular node
fixeddofs   =[1:1:2*(nely+1)]; % fixed nodes

rotM = zeros(nely,nelx); % rotation matrix

%FE_elastic(nelx,nely,x,penal,rotM);%,F,fixeddofs)
E1 = 10; % Young's mod
E2 = 5;
v1 = 0.3; % Piossons ratio
v2 = 0.25; % Piossons ratio

avgE = (E1+E2)/2;
avgV = (v1+v2)/2;
G = avgE/(2*(1+avgV));

FEA_elastic_V2(nelx,nely,x,penal,rotM,E1,E2,v1,v2,G);