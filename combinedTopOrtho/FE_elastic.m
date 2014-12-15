function [T]  = FE_elastic(nelx,nely,x,penal)



u0 =0; % value at essentail boundaries
nn = (nelx+1)*(nely+1); % number of nodes
ne = nelx*nely; % number of elements
ndof = nn*2; % Number of degrees of freedome. 2 per node. 

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
F = zeros(ndof,1);
K = zeros(ndof,ndof);
% Essential   = [1:numNodesInRow ... % bottom row
%             numNodesInRow*(numNodesInColumn-1):numNodesInRow*numNodesInColumn ... % top row
%             1:numNodesInRow :numNodesInRow*(numNodesInColumn-1) ... % left side
%             numNodesInRow:numNodesInRow:numNodesInRow*numNodesInColumn ];  % right side
        
% % Just the left side
% Essential=   1:numNodesInRow :numNodesInRow*(numNodesInColumn-1); % ... % left side
% 
% % Just the four corners
% Essential= [1 2 
%             numNodesInRow  numNodesInRow+1
%             (numNodesInColumn-1)*numNodesInRow+1     (numNodesInColumn-1)*numNodesInRow+2 
%             numNodesInColumn*numNodesInRow*2         numNodesInColumn*numNodesInRow*2-1];
        
% Essential   =[ 1:numNodesInRow*2 :2*numNodesInRow*(numNodesInColumn-1)+1 % X direction on left side is held fixed
%              1+1:numNodesInRow*2 :2*numNodesInRow*(numNodesInColumn-1)+1+1  ];  % Y direction on left side is held fixed
         
 Essential   =union( 1:2*numNodesInRow :2*(numNodesInColumn-1)*numNodesInRow+1,numNodesInRow*2) ;% X direction on left side is held fixed
               % Y direction on left side is held fixed
Essential = unique(Essential);
alldofs     = [1:ndof];
Free    = setdiff(alldofs,Essential);

% Add source node in the middle. Keep the temperature constant
%xMiddle  = floor(2*numNodesInRow/2);
%yMiddle = floor(2*numNodesInColumn/2);
%nodeNumber = 2*(numNodesInRow)*floor( (numNodesInColumn)/2);
% F(nodeNumber,1) = -1; % force at particular node

F(2*(numNodesInColumn-1)*numNodesInRow+2,1) = -10; % y direction, down


ke = elK_elastic();
  
xLoc = 1;
yLoc = 1;
% % loop over the elements
for e = 1:ne
     % compute local element vector fe
     
      % loop over local node numbers to get their node global node numbers
      for j = 1:4
          % Get the node number
          coordNodeNumber = IEN(e,j);
           % get the global X,Y position of each node and put in array
           coord(j,:) = globalPosition(coordNodeNumber,:);
      end
      
      % ----------------------
      % Pre calculated the ke matrix, so we don't have to calculated it
      % each time. 
      % ----------------------
      
      
%       % Loop over the Guassian quadrature points. We will use 2 point
%       % quadrature for a total of 4 points
%       etaRow(1,:) = [1/sqrt(3) 1/sqrt(3) -1/sqrt(3) -1/sqrt(3)];
%       zetaRow(1,:) = [1/sqrt(3) -1/sqrt(3) -1/sqrt(3) 1/sqrt(3)];
%       weight = [ 1 1 1 1];
%       
%       ke = zeros(8,8);
% %      ftemp = zeros(4,4);
%       btemp = zeros(2,4);
%      
%       for gu = 1:4
%           eta = etaRow(gu);
%           Zeta = zetaRow(gu);
%           wght = weight(gu);
%           
%           % Calculate the Shape functions.
%           N(1) = 1/4*(1-eta)*(1-Zeta);
%           N(2) = 1/4*(1+eta)*(1-Zeta);
%           N(3) = 1/4*(1+eta)*(1+Zeta);
%           N(4) = 1/4*(1-eta)*(1+Zeta);
%           
%           % B_hat (Derivative of N1 with respect to zeta and eta)
%           B_hat = 1/4*[-(1-eta) (1-eta) (1+eta) -(1+eta);
%                        -(1-Zeta) -(1+Zeta) (1+Zeta) (1-Zeta)];
%           
%           % Calculate the Jacobian
%           J=B_hat*coord;
%           
%           % Calculate the determinate
%           J_det = det(J);
%           J_transpose = transpose(J);
%           J_transpose_inv = inv(J_transpose);
%           
%           % Form the B matrix          
%          B_2by4 = J_transpose_inv*B_hat;
% 
%          % Form B, which is an 3 by 8
%          B(1,[1,3,5,7]) = B_2by4(1,1:4);
%          B(2,[2,4,6,8]) = B_2by4(2,1:4);
% 
%          B(3,[1,3,5,7]) = B_2by4(2,1:4);
%          B(3,[2,4,6,8]) = B_2by4(1,1:4);
%          
%           D = [ 1 v 0;
%               v 1 0;
%               0 0 1/2*(1-v)]*E/(1-v^2);
%           
%          
%         %  btemp = btemp +B;
%           
%           tempK = transpose(B)*D*B*J_det*wght;          
%            ke = ke + tempK;% +tempM;
%       end  
      
   %   ke

      % Insert the element stiffness matrix into the global.    
     % node = IEN(e,:);
        nodes1 = IEN(e,:);
         xNodes = nodes1*2-1;
        yNodes = nodes1*2;
    
     % NodeNumbers = union(xNodeNumbers,yNodeNumbers);
     
     % I cannot use the union, or else the order get messed up. The order
     % is important. Same in the actual topology code when you are
     % calculating the objectiv
      NodeNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
      
      % for the x location
      % The first number is the row - "y value"
      % The second number is the column "x value"
       K(NodeNumbers,NodeNumbers) = K(NodeNumbers,NodeNumbers) + x(yLoc,xLoc)^penal*ke;
       
       xLoc = xLoc+1;
       if(xLoc>nelx)
           xLoc = 1;
           yLoc = yLoc+1;
       end

end
  
K_ff = K(Free,Free);
K_fe = K(Free,Essential);
F_f = F(Free);
  
T(Free) = K_ff \ F_f;
T(Essential) = u0;
  
% disp('The temperature at each node is');
% T_column = [transpose(T),transpose(1:nn)]  
%  
% % Calcualate the heat flux
% disp('The heat flux in each element is');
%  
% qstored = zeros(ne,2);
% qMag_stored = zeros(nn,1);
% elemcenterLocations = zeros(ne,2); 
% 
% subplot(2,2,1)
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

TcontourMatrixX = zeros(numNodesInRow,numNodesInColumn);
TcontourMatrixY = zeros(numNodesInRow,numNodesInColumn);

 
  
for j = 1:numNodesInColumn % y
      rowMultiplier = j-1;
     for i = 1:numNodesInRow % x
         nodeNumber = i+numNodesInRow*rowMultiplier;
         % combined the X and Y into 1 magnitude of a vector
         node1 = nodeNumber*2-1; % minus 1, because we want to start at 1 and not 2
         node2 = node1+1;
         TcontourMatrixX(i,j) = T(node1);
         TcontourMatrixY(i,j) = T(node2);
         
         temp = sqrt(T(node1)^2+T(node2)^2);
         
         
         TcontourMatrix(i,j) = temp;
     
     end
end

if(1==0)
     subplot(1,2,1)
     % plot the coutour graph
     contour(XLocations,YLocations,TcontourMatrix);
     tti= strcat('Heat contours. Number of elements=', int2str(ne));
     title(tti);
     
    % 
    % plot the surf graph
    subplot(1,2,2)
    surf(XLocations,YLocations,TcontourMatrix);
end
 subplot(2,2,1);
  surf(XLocations,YLocations,TcontourMatrix);
  
if(1==1)
     subplot(2,2,2);
     surf(XLocations,YLocations,TcontourMatrixX);
      subplot(2,2,3);
     surf(XLocations,YLocations,TcontourMatrixY);
end
% tti= strcat('Heat surface.  Number of elements =', int2str(ne));
% title(tti);
% 
% diary off
T = transpose(T);
