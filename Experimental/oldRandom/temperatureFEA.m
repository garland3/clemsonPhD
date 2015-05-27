function temperatureFEA()

% Number  counterclock wise.
clear;
clc;
close all;
% diary('cmdwindow.txt')
% diary on
kmaterial = 1; % W/degreeC
u0 =20; % temperature at the essential boundaries


% Let the material be copper
kconductance = 401 ; % W/m  K
k = 1/kconductance;

nelx = 40;
nely = 30;


% function [U]=FE(nelx,nely,x,penal,F,fixeddofs)
[KE] = elementKtemperature(); 




K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
%F = sparse(2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
for elx = 1:nelx
  for ely = 1:nely
    n1 = (nely+1)*(elx-1)+ely; 
    n2 = (nely+1)* elx   +ely;
    edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;
  end
end
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
% F(2,1) = -1;
% fixeddofs   = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)])
alldofs     = [1:2*(nely+1)*(nelx+1)];
freedofs    = setdiff(alldofs,fixeddofs);
% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);      
U(fixeddofs,:)= 0;

% 
% elementsPerRow = 8;
% nodesInRow = elementsPerRow+2;
% 
% ne = elementsPerRow^2; % number of elements
% nn = (elementsPerRow+1)^2; % number of nodes 
% 
% count = 1;
% for j= 0:(elementsPerRow-1)       
% 
%     step=j*elementsPerRow;
%     step2=j*(elementsPerRow+1);
%     IEN(1+step,:) = [1+step2 2+step2 1+nodesInRow+step2 0+nodesInRow+step2];
%     IEN(2+step,:) = [2+step2 3+step2 2+nodesInRow+step2 1+nodesInRow+step2];
%     IEN(3+step,:) = [3+step2 4+step2 3+nodesInRow+step2 2+nodesInRow+step2];
%     IEN(4+step,:) = [4+step2 5+step2 4+nodesInRow+step2 3+nodesInRow+step2];
%     IEN(5+step,:) = [5+step2 6+step2 5+nodesInRow+step2 4+nodesInRow+step2];
%     IEN(6+step,:) = [6+step2 7+step2 6+nodesInRow+step2 5+nodesInRow+step2];
%     IEN(7+step,:) = [7+step2 8+step2 7+nodesInRow+step2 6+nodesInRow+step2];
%     IEN(8+step,:) = [8+step2 9+step2 8+nodesInRow+step2 7+nodesInRow+step2];
% 
% end
% 
% 
% XLocations=zeros((elementsPerRow+1),(elementsPerRow+1));
% YLocations=zeros((elementsPerRow+1),(elementsPerRow+1));
% 
% globalPosition = zeros(nn,2);
% globalPosition(1,:) = [ 0 0];
% for i = 1:nn     
%       column = mod(i-1,(elementsPerRow+1));
%       row = floor((i-1)/(elementsPerRow+1));
%       globalPosition(i,:)=[(column)/elementsPerRow (row)/elementsPerRow];
% 
%       XLocations(row+1,column+1)= globalPosition(i,1);
%       YLocations(row+1, column+1) = globalPosition(i,2);
% end
% 
% Essential   = [1:9    72:81    1:9:73    9:9:81 ]; 
% Essential = unique(Essential);
% alldofs     = [1:nn];
% Free    = setdiff(alldofs,Essential);
%     
% 
% 
% 
%  
%  % List the node numbers start-end where there is flux
% %listOfFluxNodes = [5 9;
% %                    9 12; 12 14; 14 15; 15 13; 13 10; 10 6; 6 1];
%  
%      
% F = zeros(nn,1);
% K = zeros(nn,nn);
% B_stored = zeros(ne,8);
%  
% % loop over the elements
% for e = 1:ne
function [KE]=elementKtemperature
% compute local element vector fe

 % loop over local node numbers to get their node global node numbers
 
 % Each box is 1 by 1
coord(1,:) = [0 0]
coord(2,:) = [1 0]
coord(3,:) = [1 1]
coord(4,:) = [0 1]
 
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
 ke
 %B_stored(e,1:8)= btemp(1:8);

% Insert the element stiffness matrix into the global.
 % node = IEN(e,:);

 % K(node,node) = K(node,node) +ke;

  % Now do something similar for F
 % Calculate g at the nodes of this element
%  g = sin(pi*coord(:,1)).*sin(pi*coord(:,2));  
% 
%  % multiply by the integrated ftemp
%  fe = ftemp*g;
% 
%  % Add the f_e element to the F global matrix
%  node = IEN(e,:);
%  F(node) = F(node) + fe;
     
% end
%  
%  
% K_ff = K(Free,Free);
% K_fe = K(Free,Essential);
% F_f = F(Free);
%  
% T(Free) = K_ff \ F_f;
% T(Essential) = u0;
%  
% disp('The temperature at each node is');
% T_column = [transpose(T),transpose(1:nn)]
%  
%  
% % Calcualate the heat flux
% disp('The heat flux in each element is');
%  
% qstored = zeros(ne,2);
% qMag_stored = zeros(nn,1);
% elemcenterLocations = zeros(ne,2);
% 
% 
% subplot(2,2,1)
%  
% 
% 
% % loop over the elements
% for e = 1:ne
%     
%     coord = zeros(3,2);
%     xsum = 0;
%     ysum = 0;
%     % loop over local node numbers to get their node global node numbers
%     for j = 1:4
%         % Get the node number
%         coordNodeNumber = IEN(e,j);
%          % get the global X,Y position of each node and put in array
%          coord(j,:) = globalPosition(coordNodeNumber,:);
%          local_t(j) = T(coordNodeNumber);
%          xsum = xsum+coord(j,1);
%          ysum = ysum+coord(j,2);
%     end
%     elemcenterLocations(e,:) = [xsum/4 ysum/4];
%     
%     % see version 3 of notes page 13 Also, see version 5 of the notes page
%     % 27
%    
%     eta = 0; Zeta = 0; % We are at the center, so both are zero
%    
%     % B_hat (Derivative of N1 with respect to zeta and eta)
%      B_hat = 1/4*[-(1-eta) (1-eta) (1+eta) -(1+eta);
%                   -(1-Zeta) -(1+Zeta) (1+Zeta) (1-Zeta)];
% 
%      % Calculate the Jacobian
%      J=B_hat*coord;
% 
%      % Calculate the determinate
%      %J_det = det(J);
%      J_transpose = transpose(J);
%      J_transpose_inv = inv(J_transpose);
% 
%      % Form the B matrix
%      B = J_transpose_inv*B_hat;
%     
%     qLocal = -kmaterial*B*local_t';
%     qstored(e,:) = qLocal';
%     
%     qMag = (qLocal(1)^2+qLocal(2)^2)^(1/2);
%     qMag_stored(e) =qMag;
%     
%     % plot the element outline
%     hold on
%     coord(5,:) = coord(1,:); 
%     plot(coord(:,1),coord(:,2));    
%     
% end
% 
% quiver(elemcenterLocations(:,1),elemcenterLocations(:,2),qstored(:,1),qstored(:,2))
% %xlabel('radial distance, meters') % y-axis label
% %ylabel('Height') % x-axis label
% tti= strcat('Flux from each element . Number of elements=', int2str(ne));
% title(tti);
% hold off
% 
% q_mags = [qMag_stored,transpose(1:nn)]
% 
% subplot(2,2,2)
%  TcontourMatrix = zeros((elementsPerRow+1),(elementsPerRow+1));
% 
%  
%  for j = 1:(elementsPerRow+1)
%      for i = 1:(elementsPerRow+1)
%          TcontourMatrix(i,j) = T((j-1)*(elementsPerRow+1)+i);
%      
%      end
%  end
%  
%  % plot the coutour graph
% contour(XLocations,YLocations,TcontourMatrix);
% tti= strcat('Heat contours. Number of elements=', int2str(ne));
% title(tti);
% 
% % plot the surf graph
% subplot(2,2,3)
% surf(XLocations,YLocations,TcontourMatrix);
% 
% tti= strcat('Heat surface.  Number of elements =', int2str(ne));
% title(tti);
% 
% diary off