function [KE]=elK_elastic

% Change this if statement to be not true if you need to see the
% calculations
if(1==0)
    KE =   [ 1.9556    0.6667   -1.1556   -0.1333   -0.9778   -0.6667    0.1778    0.1333;
            0.6667    1.9556    0.1333    0.1778   -0.6667   -0.9778   -0.1333   -1.1556;
           -1.1556    0.1333    1.9556   -0.6667    0.1778   -0.1333   -0.9778    0.6667;
           -0.1333    0.1778   -0.6667    1.9556    0.1333   -1.1556    0.6667   -0.9778;
           -0.9778   -0.6667    0.1778    0.1333    1.9556    0.6667   -1.1556   -0.1333;
           -0.6667   -0.9778   -0.1333   -1.1556    0.6667    1.9556    0.1333    0.1778;
            0.1778   -0.1333   -0.9778    0.6667   -1.1556    0.1333    1.9556   -0.6667;
            0.1333   -1.1556    0.6667   -0.9778   -0.1333    0.1778   -0.6667    1.9556];
    
    

else
% Do the calculations 
E1 = 10; % Young's mod
E2 = 5;
v1 = 0.3; % Piossons ratio
v2 = 0.25; % Piossons ratio

avgE = (E1+E2)/2;
avgV = (v1+v2)/2;
G = avgE/(2*(1+avgV));

D = 1/(1-v1*v2)*[E1 v2*E1 0;
                v1*E2 E2 0;
                0     0  G*(1-v1*v2)];
            
            
  D = [ 1 avgV 0;
          avgV 1 0;
          0 0 1/2*(1-avgV)]*avgE/(1-avgV^2);

% Rotation matrix
% from a book. I need to verify this myself
% Dnew = transpose(T)*D*T
% T = [cos(theta)^2 sin(theta)^2 2*sin(theta)*cos(theta);
%     sin(theta)^2 cos(theta)^2 -2*sin(theta)*cos(theta);
%     2*sin(theta)*cos(theta) -2*sin(that)*cos(theta) cos(theta)^2-sin(theta)^2]

% 1 by 1 square
coord = [0 0;
        1 0;
        1 1;
        0 1];
    
%     % try sigmund's square, down is positive, then I get the same answer. 
%     coord = [0 0;
%         1 0;
%         1 1;
%         0 1];
    
      etaRow(1,:) = [1/sqrt(3) 1/sqrt(3) -1/sqrt(3) -1/sqrt(3)];
      zetaRow(1,:) = [1/sqrt(3) -1/sqrt(3) -1/sqrt(3) 1/sqrt(3)];
      weight = [ 1 1 1 1];
      
      ke = zeros(8,8);
%      ftemp = zeros(4,4);
      btemp = zeros(2,4);
     
      % Loop over the guass points
      for gu = 1:4
          eta = etaRow(gu);
          Zeta = zetaRow(gu);
          wght = weight(gu);
          
          % Calculate the Shape functions.
          N(1) = 1/4*(1-eta)*(1-Zeta);
          N(2) = 1/4*(1+eta)*(1-Zeta);
          N(3) = 1/4*(1+eta)*(1+Zeta);
          N(4) = 1/4*(1-eta)*(1+Zeta);
          
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
      KE = ke;
end
