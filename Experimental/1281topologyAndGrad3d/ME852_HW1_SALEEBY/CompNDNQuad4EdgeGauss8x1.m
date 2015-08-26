%--------------------------------------------------------------------
% Compute N(x,y), dNdx(x,y), dNdy(x,y) at the Gauss points on the edges 
% of a square master element
%
%       --2--1--
%      |        |
%      3        8
%      |        |
%      4        7
%      |        |
%       --5--6--
% 
%--------------------------------------------------------------------
function [N, Nx, Ny, Nz]=CompNDNQuad4EdgeGauss8x1()

  [gauss_points, gauss_weights]=GetQuadEdgeGauss8x1();
  
  N=zeros(8,24);
  Nx=zeros(8,24);
  Ny=zeros(8,24);
  Nz=zeros(8,24);
  
  for j=1:24
    x=gauss_points(j,1);
    y=gauss_points(j,2);
    z=gauss_points(j,3);
    [Np,Npx,Npy,Npz]=CompNDNatPointQuad4(x,y,z);
    N(:,j)=Np;
    Nx(:,j)=Npx;
    Ny(:,j)=Npy;
    Nz(:,j)=Npz;
  end
