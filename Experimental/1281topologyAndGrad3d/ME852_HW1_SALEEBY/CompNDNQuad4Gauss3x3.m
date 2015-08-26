
function [N,Nx,Ny,Nz]=CompNDNQuad4Gauss3x3()
  
  [gauss_points, gauss_weights]=GetQuadGauss3x3();
  
  N=zeros(8,8);
  Nx=zeros(8,8);
  Ny=zeros(8,8);
  Nz=zeros(8,8);

  
  for j=1:8
    x=gauss_points(j,1);
    y=gauss_points(j,2);
    z=gauss_points(j,3);
    [Np,Npx,Npy,Npz]=CompNDNatPointQuad4(x,y,z);
    N(:,j)=Np;
    Nx(:,j)=Npx;
    Ny(:,j)=Npy;
    Nz(:,j)=Npz;
  end
