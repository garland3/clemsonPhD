function [KE]=elementKtemperature()
kmaterial = 1;
% compute local element vector fe

 % loop over local node numbers to get their node global node numbers
 
 % Each box is 1 by 1
 
% the top left corner is 1
% bttom left is 2
% bottom right is 3
% top right is 4
% 
% down is positive y
% right is positive x
coord(1,:) = [0 0]
coord(2,:) = [0 1]
coord(3,:) = [1 1]
coord(4,:) = [1 0]
 
%  
%  for j = 1:4
%      % Get the node number
%      coordNodeNumber = IEN(e,j);
%       % get the global X,Y position of each node and put in array
%       coord(j,:) = globalPosition(coordNodeNumber,:);
%  end

 % See notes page 18 Calculate the Jacobian
 x1 = coord(1,1);
 x2 = coord(2,1);
 x3 = coord(3,1);
 x4 = coord(4,1);

 y1 = coord(1,2);
 y2 = coord(2,2);
 y3 = coord(3,2);
 y4 = coord(4,2);


 % Loop over the Guassian quadrature points. We will use 2 point
 % quadrature for a total of 4 points
 etaRow(1,:) = [1/sqrt(3) 1/sqrt(3) -1/sqrt(3) -1/sqrt(3)];
 zetaRow(1,:) = [1/sqrt(3) -1/sqrt(3) -1/sqrt(3) 1/sqrt(3)];
 weight = [ 1 1 1 1];

 ke = zeros(4,4);
 ftemp = zeros(4,4);
 btemp = zeros(2,4);

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
     J_transpose = transpose(J);
     J_transpose_inv = inv(J_transpose);

     % Form the B matrix
     B = J_transpose_inv*B_hat;
     btemp = btemp +B;

     tempK = transpose(B)*kmaterial*B*J_det*wght;

     % now calculate the m term
     % tempM = transpose(N)*10*N*J_det*wght; % 


      ke = ke + tempK;%+tempM;

     % ----------------------------- Now calculate the f vector

     % Calculate the x and y position of this guass point
    % x = N(1)*x1+N(2)*x2+N(3)*x3+N(4)*x4; y =
    % N(1)*y1+N(2)*y2+N(3)*y3+N(4)*y4;

     % calculate the g points
    % g(1) = sin(pi*x)*sin(pi*y);
     temp = transpose(N)*2*pi^2*N*J_det*wght; % *g??
     ftemp = ftemp+ temp;              

 end  
 KE=ke