nelx = 20;
nely = 10;
x = ones(nely,nelx);
F = sparse(2*(nely+1)*(nelx+1),1);
penal=1;


nodeNumber = 2*(nely+1)*floor( (nelx+1)/2);
F(nodeNumber,1) = -10; % force at particular node
fixeddofs   =[1:1:2*(nely+1)]; % fixed nodes

FE_elastic(nelx,nely,x,penal);%,F,fixeddofs)