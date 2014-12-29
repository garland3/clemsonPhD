function [ ke ] = ke_elementV2( theta,E1,E2,v1,v2,G )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% 1 by 1 square
coord = [0 0;
        1 0;
        1 1;
        0 1];

    




            
Dtemp = 1/(1-v1*v2)*[E1 v2*E1 0;
                v1*E2 E2 0;
                0     0  G*(1-v1*v2)];
            
% Dnew = transpose(T)*D*T
% See pdf Elasticity_of_Composite_materials

TransFormMatrix =  [cos(theta)^2 sin(theta)^2 2*sin(theta)*cos(theta);
                     sin(theta)^2 cos(theta)^2 -2*sin(theta)*cos(theta);
                     -2*sin(theta)*cos(theta) 2*sin(theta)*cos(theta) cos(theta)^2-sin(theta)^2]; % page 57, equation 2.32
            
D  = transpose(TransFormMatrix)*Dtemp*TransFormMatrix;      % page 59 equattion 2.34                
% ----------------------
% Calculate the element stiffness matrix
% each time. 
% ----------------------
  etaRow(1,:) = [1/sqrt(3) 1/sqrt(3) -1/sqrt(3) -1/sqrt(3)];
  zetaRow(1,:) = [1/sqrt(3) -1/sqrt(3) -1/sqrt(3) 1/sqrt(3)];
  weight = [ 1 1 1 1];

  ke = zeros(8,8);


  % Loop over the guass points
  for gu = 1:4
      eta = etaRow(gu);
      Zeta = zetaRow(gu);
      wght = weight(gu);

      % page 9 of solutions for hw 8 as a reference
      % B_hat (Derivative of N1 with respect to zeta and eta)
      B_hat = 1/4*[-(1-eta) (1-eta) (1+eta) -(1+eta);
                   -(1-Zeta) -(1+Zeta) (1+Zeta) (1-Zeta)];

      % Calculate the Jacobian
      J=B_hat*coord;

      % Calculate the determinate
      J_det = det(J);        
      % J_inverse = inv(J);         

      % Form the B matrix          
      % B_2by4 = J_inverse*B_hat;
      B_2by4_v2 = J\B_hat;

     % Form B, which is an 3 by 8
     B = zeros(3,8);
     B(1,[1,3,5,7]) = B_2by4_v2(1,1:4);
     B(2,[2,4,6,8]) = B_2by4_v2(2,1:4);

     B(3,[1,3,5,7]) = B_2by4_v2(2,1:4);
     B(3,[2,4,6,8]) = B_2by4_v2(1,1:4); 

    

     tempK = transpose(B)*D*B*J_det*wght;          
     ke = ke + tempK; 
  end  
      
