
%--------------------------------------------------------------------

function B= CompB4x8Quad4Gauss2x2(g, Nx, Ny, Nz)

  B=zeros(9,24);
  
  for i=1:8
    B(1,3*i-2)=Nx(i,g);
    B(4,3*i-1)= Nx(i,g);
    B(7,3*i)= Nx(i,g);
    
    B(2,3*i-2)=Ny(i,g);
    B(5,3*i-1)= Ny(i,g);
    B(8,3*i)= Ny(i,g);

    B(3,3*i-2)=Nz(i,g);
    B(6,3*i-1)= Nz(i,g);
    B(9,3*i)= Nz(i,g);

  end
