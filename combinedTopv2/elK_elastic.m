function [KE]=elK_elastic

% Change this if statement to be not true if you need to see the
% calculations
if(1==1)
    KE =    1.0e+07 * [1.4305    0.5107   -0.8683   -0.0515   -0.7152   -0.5107    0.1531    0.0515;
                        0.5107    1.4305    0.0515    0.1531   -0.5107   -0.7152   -0.0515   -0.8683;
                       -0.8683    0.0515    1.4305   -0.5107    0.1531   -0.0515   -0.7152    0.5107;
                       -0.0515    0.1531   -0.5107    1.4305    0.0515   -0.8683    0.5107   -0.7152;
                       -0.7152   -0.5107    0.1531    0.0515    1.4305    0.5107   -0.8683   -0.0515;
                       -0.5107   -0.7152   -0.0515   -0.8683    0.5107    1.4305    0.0515    0.1531;
                        0.1531   -0.0515   -0.7152    0.5107   -0.8683    0.0515    1.4305   -0.5107;
                        0.0515   -0.8683    0.5107   -0.7152   -0.0515    0.1531   -0.5107    1.4305];
    
    

else
    % Do the calculations 
    E = 29007547; % psi,  Young's mod
    v = 0.29; % Piossons ratio
    % G = E/(2*(1+v));

    D = [ 1 v 0;
          v 1 0;
          0 0 1/2*(1-v)]*E/(1-v^2);


    % 1 by 1 square
    coord = [0 0;
            1 0;
            1 1;
            0 1];

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
  KE = ke;
end
