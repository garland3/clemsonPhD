%--------------------------------------------------------------------
% Compute material matrix (plane stress)
%--------------------------------------------------------------------

function C= CompCPlaneStress(materials)

  C=zeros(6,6);
  E=materials(1);  % Young's modulus
  nu=materials(2); % Poisson's ratio

  t = E/((1+nu)*(1-2*nu));
  C(1,1)=(1-nu);
  C(2,2)=(1-nu);
  C(3,3)=(1-nu);
  
  C(1,2)=nu;
  C(1,3)=nu;
  C(2,1)=nu;
  C(2,3)=nu;
  C(3,1)=nu;
  C(3,2)=nu;

  C(4,4)=(1-2*nu)/2;
  C(5,5)=(1-2*nu)/2;
  C(6,6)=(1-2*nu)/2;
  C=t*C;
end