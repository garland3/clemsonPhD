function [KE, kExpansion_bar, B_total, F_meso]=elK_elastic(E,v, G, strain,Dgiven)

% Change this if statement to be not true if you need to see the
% calculations
if(1==1)
  
    % Do the calculations 
    %E = 1; % Young's mod


    % E = materialProp.E_material1;
    % v = materialProp.v; % Piossons ratio
    % G = materialProp.G;
    
    % Plane stress problem. 
      if(isempty(Dgiven))
        D = [ 1 v 0;
            v 1 0;
            0 0 1/2*(1-v)]*E/(1-v^2);
      else
          D=Dgiven;
      end


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

     % etaRow(1,:) = [1/sqrt(3) 1/sqrt(3) -1/sqrt(3) -1/sqrt(3)];
     etaRow(1,:)= [0.5774    0.5774   -0.5774   -0.5774];
     % zetaRow(1,:) = [1/sqrt(3) -1/sqrt(3) -1/sqrt(3) 1/sqrt(3)];
     zetaRow(1,:) =  [  0.5774   -0.5774   -0.5774    0.5774];
     weight = [ 1 1 1 1];

     ke = zeros(8,8);
     F_meso = zeros(8,1);
     kExpansion_bar = 0;
%      ftemp = zeros(8,1);
    %  btemp = zeros(2,4);
     B_total = zeros(3,8);

      % Loop over the guass points
      for gu = 1:4
          eta = etaRow(gu);
          Zeta = zetaRow(gu);
          wght = weight(gu);

%           % Calculate the Shape functions.
%           N(1) = 1/4*(1-eta)*(1-Zeta);
%           N(2) = 1/4*(1+eta)*(1-Zeta);
%           N(3) = 1/4*(1+eta)*(1+Zeta);
%           N(4) = 1/4*(1-eta)*(1+Zeta);

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
         B_total = B_total+ B*J_det*wght;
         

         tempK = transpose(B)*D*B*J_det*wght;          
         ke = ke + tempK; 
         
         
       
         
         % Calcualte the thermal expansion matrix
         f_bar_temp =  transpose(B)*D*transpose([1 1  0])*J_det*wght;
         kExpansion_bar = kExpansion_bar+f_bar_temp;
      end  
      
      KE = ke;
      
    % calculate the load force for a strain 
     if(~isempty(strain))
      F_meso = transpose(B_total)*D*transpose(strain);
     end
      
     
      

elseif(1==2)      
    KE =   [ 1.9556    0.6667   -1.1556   -0.1333   -0.9778   -0.6667    0.1778    0.1333;
            0.6667    1.9556    0.1333    0.1778   -0.6667   -0.9778   -0.1333   -1.1556;
           -1.1556    0.1333    1.9556   -0.6667    0.1778   -0.1333   -0.9778    0.6667;
           -0.1333    0.1778   -0.6667    1.9556    0.1333   -1.1556    0.6667   -0.9778;
           -0.9778   -0.6667    0.1778    0.1333    1.9556    0.6667   -1.1556   -0.1333;
           -0.6667   -0.9778   -0.1333   -1.1556    0.6667    1.9556    0.1333    0.1778;
            0.1778   -0.1333   -0.9778    0.6667   -1.1556    0.1333    1.9556   -0.6667;
            0.1333   -1.1556    0.6667   -0.9778   -0.1333    0.1778   -0.6667    1.9556];
else
    
    


end
