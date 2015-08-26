
function [N,Nx,Ny,Nz]=CompNDNatPointQuad4(x,y,z)
  
  N=zeros(8,1);
  Nx=zeros(8,1);
  Ny=zeros(8,1);
  Nz=zeros(8,1);

%   master_nodes=[1 1; -1 1; -1 -1; 1 -1];
% master_nodes=[-1 -1 -1; 1 -1 -1; -1 1 -1; 1 1 -1;
%-1 -1 1; 1 -1 1; -1 1 1; 1 1 1];  

%   master_nodes=[1 1 1; -1 1 1; -1 -1 1; 1 -1 1;
%       1 1 -1; -1 1 -1; -1 -1 -1; 1 -1 -1];
  master_nodes=[1 1 -1; -1 1 -1; -1 -1 -1; 1 -1 -1;
      1 1 1; -1 1 1; -1 -1 1; 1 -1 1];
  
  for i=1:8
    nx=master_nodes(i,1); 
    ny=master_nodes(i,2);
    nz=master_nodes(i,3);
    N(i)=(1.0 + nx*x)*(1.0 + ny*y)*(1.0+nz*z)/8.0;
    Nx(i)= nx*(1.0 + ny*y)*(1.0 + nz*z)/8.0;
    Ny(i)= ny*(1.0 + nx*x)*(1.0 + nz*z)/8.0;
    Nz(i)= nz*(1.0 + nx*x)*(1.0 + ny*y)/8.0;
  end
