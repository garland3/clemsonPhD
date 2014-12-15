function top()
sizeX =30; % number of elements in X direction
sizeY = 15; % number of elements in Y direction
TargetVolFraction = 0.3;
Penalty = 1;
rmin = 2;

% boundary conditions for fea
nelx = sizeX
nely = sizeY;
F = sparse(2*(nely+1)*(nelx+1),1);
E1 = 10;
E2 = 1;


% Default loading
% F(2,1) = -1; % force at particular node
% fixeddofs   = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]); % fixed nodes

% Middle Bottom load, left side is fixed
nodeNumber = 2*(nely+1)*floor( (nelx+1)/2)
F(nodeNumber,1) = -10; % force at particular node
fixeddofs   =[1:1:2*(nely+1)]; % fixed nodes

%ActualtopTemperature(sizeX,sizeY,TargetVolFraction,Penalty,rmin,F,fixeddofs,E1,E2)   

Actualtop(sizeX,sizeY,TargetVolFraction,Penalty,rmin,F,fixeddofs)

