function [T]  = temperatureFEA2(nelx,nely,x,penal, doPlot, iterationsPerPlot, iteration)

kmaterial = 1; % W/degreeC
u0 =0; % value at essentail boundaries
nn = (nelx+1)*(nely+1); % number of nodes
ne = nelx*nely; % number of elements

% ------------------------------------
% -- Make a mapping between elements and global matrixes
% ------------------------------------
% IEN holds the node numbers for each element. 
% Each row is a new element
% The first column is element 1's global node number
% Second column is elemnt 2's global node number
%  and ....

count = 1;
elementsInRow = nelx+1;
IEN = zeros(nn,4);
% Each row, so nely # of row
for i = 1:nely
     rowMultiplier = i-1;
    % Each column, so nelx # of row
    for j= 1:nelx        
        IEN(count,:)=[rowMultiplier*elementsInRow+j,
                      rowMultiplier*elementsInRow+j+1,
                      (rowMultiplier +1)*elementsInRow+j+1,
                       (rowMultiplier +1)*elementsInRow+j];
        count = count+1;
    end
end

numNodesInRow = nelx+1;
numNodesInColumn = nely+1;
XLocations=zeros(numNodesInRow,numNodesInColumn);
YLocations=zeros(numNodesInRow,numNodesInColumn);

% Find and store the global positions of each node
% Each element square is 1 by 1 units
% Store both the X and Y positions
globalPosition = zeros(nn,2);
count = 1;
for i = 1:numNodesInColumn  % y
    for j= 1:numNodesInRow % x
        globalPosition(count,:) = [j-1 i-1];
        count = count+1;        
        XLocations(j,i) = j-1;
        YLocations(j,i) = i-1;
    end
end 

% Specifiy the constrained nodes where there are essential boundary
% conditions
F = zeros(nn,1);
K = zeros(nn,nn);

row = nelx+1;
column = nely+1;
% Just the left side in the middle
quartY = ceil(nely/4);
Essential=   [1 row (column-1)*row+1 column*row] ; % the 4 corners

Essential = unique(Essential);
alldofs     = [1:nn];
Free    = setdiff(alldofs,Essential);
%F = ones(nn,1); % add a constant heat source everywhere
F([ceil(row/2)+(ceil(column/2)*row) (ceil(row/2)+1)+(ceil(column/2)*row)]) =  20; % heat source in the middle

  
xLoc = 1;
yLoc = 1;
% loop over the elements
for e = 1:ne
    
      % Get the precalculated element stiffness matrix. 
      ke = elementK_heat();
      
      
      % Insert the element stiffness matrix into the global.    
      node = IEN(e,:);
      
      % for the x location
      % The first number is the row - "y value"
      % The second number is the column "x value"
       K(node,node) = K(node,node) + x(yLoc,xLoc)^penal*ke;
       
       xLoc = xLoc+1;
       if(xLoc>nelx)
           xLoc = 1;
           yLoc = yLoc+1;
       end

end
  
K_ff = K(Free,Free);
% K_fe = K(Free,Essential);
F_f = F(Free);
  
T(Free) = K_ff \ F_f;
T(Essential) = u0;

  
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

TcontourMatrix = zeros(numNodesInRow,numNodesInColumn);

 
  
for j = 1:numNodesInColumn % y
      rowMultiplier = j-1;
     for i = 1:numNodesInRow % x
         nodeNumber = i+numNodesInRow*rowMultiplier;
         TcontourMatrix(i,j) = T(nodeNumber);
     
     end
end

if(doPlot==1 && mod(iteration,iterationsPerPlot) ==0)
     subplot(1,2,1)
     % plot the coutour graph
     contour(XLocations,YLocations,TcontourMatrix);
     tti= strcat('Heat contours. Number of elements=', int2str(ne));
     title(tti);
     
    % 
    % plot the surf graph
%     subplot(1,2,2)
%     surf(XLocations,YLocations,TcontourMatrix);
end
% subplot(1,2,1);
%  surf(XLocations,YLocations,TcontourMatrix);
% tti= strcat('Heat surface.  Number of elements =', int2str(ne));
% title(tti);
% 
% diary off
T = transpose(T);
