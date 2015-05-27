% FEA
% Anthony G and Marus Yonder
% 
% Take home exam
% Problem #1 and 2

% clear; clc;  close all
problem = 1
;
a = 4; % in
d = 7.5; % in
L = 10; % in
t = 0.1; % in 
h = 1; % in

FappliedLoad = 0.25; % lb. 

% When plotting the diplaced elements, exaggerate the displacements by
% this scale
multiplierScale = 1000;

% Assume AISI 1020 Steel
% Material properties are from Solidworks material repository
E = 29007547; % psi,  Young's mod
v = 0.29; % Piossons ratio
% G = E/(2*(1+v));

% The density of the steel is given in the problem
density = 0.28; % lb/in^3

% Handle the problem 1 mapping between nodes and elements and node
% locations
if(problem ==1)
    % ==========================================================
    % ----------------------------------------------------------
    % Problem 1
    % ----------------------------------------------------------
    % ==========================================================
    
    
    nelx = 100; % Number of elements in the x direction
    nely = 10; % Number of elements in the y direction

  
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
    % 
    % Element nodes are as follows
    %
    %  4 ---- 3
    %  |      |
    %  |      |
    %  1 ---- 2
    %
    %
    % Let the x direction be the first dof and y direction the second
    %  
    %   y = 2
    %   /\
    %   |
    %   |
    %   *----> x = 1
    %
    %
    %
    % ---------------------------------------------
    % Global matrix nodes are as follows
    % row = nelx+1
    % col = nely+1
    %
    % col*row+1-col*row+2-col*row+3-col*row+4...col*row+row-1-col*row+row
    % .            .         .         .             .            .
    % .            .         .         .             .            .
    % .            .         .         .             .            .
    % |            |         |         |             |            |
    % 2*row+1 - 2*row+2 - 2*row+3 - 2*row+4 ... 2*row+row-1 - 2*row+row
    % |            |         |         |             |            |
    % 1*row+1 - 1*row+2 - 1*row+3 - 1*row+4 ... 1*row+row-1 - 1*row+row
    % |            |         |         |             |            |
    % 0*row+1 - 0*row+2 - 0*row+3 - 0*row+4 ... 0*row+row-1 - 0*row+row


    count = 1;
    numNodesInRow = nelx+1;
    numNodesInColumn = nely+1;
    IEN = zeros(nn,4); % Index of element nodes (IEN)
    % Each row, so nely # of row
    for i = 1:nely
         rowMultiplier = i-1;
        % Each column, so nelx # of row
        for j= 1:nelx        
            IEN(count,:)=[rowMultiplier*numNodesInRow+j, ...
                          rowMultiplier*numNodesInRow+j+1, ...
                          (rowMultiplier +1)*numNodesInRow+j+1, ...
                           (rowMultiplier +1)*numNodesInRow+j];
            count = count+1;
        end
    end

    % Find and store the global positions of each node
    % Each element rectangle. The aspect ratio is determined by the 
    % number off elements in the X and Y directions and the height and width of
    % the beam
    %
    % Store both the X and Y positions
    globalPosition = zeros(nn,2);
    count = 1;
    for i = 1:numNodesInColumn  % y
        for j= 1:numNodesInRow % x
            xloc = (j-1)*L/(numNodesInRow-1);
            yloc = (i-1)*h/(numNodesInColumn-1);

            globalPosition(count,:) = [xloc yloc];
            count = count+1; 
        end
    end 
else
    %% 
    
    % ==========================================================
    % ----------------------------------------------------------
    % Problem 2
    % ----------------------------------------------------------
    % ==========================================================
    refineMesh = 0; % set equal to 1 for a more detailed mesh
  
    % Get the node to element mapping from the Abaqus file
    IEN = readElementCordinateMap(refineMesh) ;      
    
    % Get the node locations from the Abaqus file
    globalPosition = readNodeCordinates(refineMesh) ;     
    
     nn = size(globalPosition,1);
     ne = size(IEN,1);
end
   
% ---------------------------------------------
% Initialize the global force, k, and B_stored matrixes
%
% ---------------------------------------------
ndof = nn*2; % Number of degrees of freedome. 2 per node. 
F2 = zeros(ndof,1);
K = zeros(ndof,ndof);
 
% ---------------------------------------------
% Find the essential boundary conditions
%
% Find the top right node and apply the point force
% ---------------------------------------------
%
% Search for the nodes on the far left to set as essential boundary
% conditions

u0 = 0; % value at essential boundaries


for i = 1:nn
    % get the xy locations as a row
    xy = globalPosition(i,:);
    xLoc = xy(1); yLoc = xy(2);
    
    if( xLoc ==0)
        % exists returns 1 if Essential2 is a varriable
        temp = exist('Essential2', 'var');
        
        if(temp==1)
             Essential2 = [Essential2 2*i 2*i-1]; %#ok<AGROW>
        else
            Essential2 =[2*i-1 2*i];
        end
    end
    
    % if in the top right corner, add the point force
    if(yLoc == 1 && xLoc == 10)
         F2(i*2) = -FappliedLoad;
    end
        
end
Essential2 = unique(Essential2);
% ---------------------------------------------
% Generate the Free dof matrix
% ---------------------------------------------
alldofs     = 1:ndof;
Free    = setdiff(alldofs,Essential2);
    

% ---------------------------------------------
%     Generate the local k and f matrixes
%     Add them to the global matrix
% ---------------------------------------------

% plane stress
% Define D outside the loop
  D = [ 1 v 0;
      v 1 0;
      0 0 1/2*(1-v)]*E/(1-v^2);
  
% % loop over the elements
for e = 1:ne         
      % loop over local node numbers to get their node global node numbers
      coord = zeros(4,2);
      for j = 1:4
          % Get the node number
          nodeNumber = IEN(e,j);
           % get the global X,Y position of each node and put in array
           coord(j,:) = globalPosition(nodeNumber,:);
      end
      
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
      
      % Calculate the body force term
      % 1. Find the area of the element
      % - evalute the jacobian at zeta = eta = 0, 
      % - then det(J)*4 = Area
      % 2. Find the total force on the element by extruding the area and
      % multiplying by the density
      % 3. Equally apply the load to each node in the negative Y dof
      
       % B_hat (Derivative of N1 with respect to zeta and eta)
       Zeta = 0; eta = 0;
       B_hat = 1/4*[-(1-eta) (1-eta) (1+eta) -(1+eta);
                       -(1-Zeta) -(1+Zeta) (1+Zeta) (1-Zeta)];
          
       % Calculate the Jacobian
       J=B_hat*coord;
        
       % Calculate the Area as determinate(Jacobian)*4
       area =  det(J)*4;
       volume = area*t; % t is the thickness
       bodyForce = volume*density;
       
       
       f_element = zeros(8,1);
       f_element(2:2:8,1) = -bodyForce/4;
      
     
     % Insert the element stiffness matrix into the global.    
     nodes1 = IEN(e,:);
     xDOF = nodes1*2-1;
     yDOF = nodes1*2;
    
     % I cannot use the union, or else the order get messed up. The order
     % is important. Same in the actual topology code when you are
     % calculating the objectiv
      dofNumbers = [xDOF(1) yDOF(1) xDOF(2) yDOF(2) xDOF(3) yDOF(3) xDOF(4) yDOF(4)];
      
      % multiply t*ke and add to global matrix. t = thickness or t_z
      K(dofNumbers,dofNumbers) = K(dofNumbers,dofNumbers) + t*ke;
      F2(dofNumbers,1) = F2(dofNumbers,1) +f_element;        
end
  
K_ff = K(Free,Free);
% K_fe = K(Free,Essential2);
F_f = F2(Free);
  
U(Free) = K_ff \ F_f;
U(Essential2) = u0;  

stress_stored = zeros(ne,3);

XYStressLocationsStored = zeros(ne,1);
vonM_stored = zeros(ne,1);

elemcenterLocations = zeros(ne,2);

% loop over the elements
 %subplot(2,2,1)
 for e = 1:ne
    arrayCoordNumber = zeros(1,4);
     xsum = 0;
     ysum = 0;
     % loop over local node numbers to get their node global node numbers
     local_u = zeros(4,1);
     for j = 1:4
         % Get the node number
         nodeNumber = IEN(e,j);
         arrayCoordNumber(j) = nodeNumber;
          % get the global X,Y position of each node and put in array
          coord(j,:) = globalPosition(nodeNumber,:);
          local_u(j,:) = U(nodeNumber);
          xsum = xsum+coord(j,1);
          ysum = ysum+coord(j,2);
     end
     
     elemcenterLocations(e,:) = [xsum/4 ysum/4];     
         
     % see version 3 of notes page 13 Also, see version 5 of the notes page
     % 27
     % Calculate stress and strains in the center of each elements
     eta = 0; Zeta = 0; % We are at the center, so both are zero

     % B_hat (Derivative of N1 with respect to zeta and eta)
      B_hat = 1/4*[-(1-eta) (1-eta) (1+eta) -(1+eta);
                   -(1-Zeta) -(1+Zeta) (1+Zeta) (1-Zeta)];

      % Calculate the Jacobian
      J=B_hat*coord;     

      % Form the B matrix          
      %B_2by4 = J_inverse*B_hat;
      B_2by4_v2 = J\B_hat;

     % Form B, which is an 3 by 8
      B = zeros(3,8);
      B(1,[1,3,5,7]) = B_2by4_v2(1,1:4);
      B(2,[2,4,6,8]) = B_2by4_v2(2,1:4);

      B(3,[1,3,5,7]) = B_2by4_v2(2,1:4);
      B(3,[2,4,6,8]) = B_2by4_v2(1,1:4);
     
     % Get the node numbers, and then the dofs
     nodes1 = IEN(e,:);
     xDOF = nodes1*2-1;
     yDOF = nodes1*2;
     dofNumbers = [xDOF(1) yDOF(1) xDOF(2) yDOF(2) xDOF(3) yDOF(3) xDOF(4) yDOF(4)];
     
     d_local = U(dofNumbers);
     
     % Strain is B*d, NOTE: d is transposed
     strain = B*d_local';    
     stress = D*strain;
     
     % Store the transpose, to make things work nice. 
     stress_stored(e,:) = stress';
     
     %avg = (stress(1)+stress(2))/2;
     %R = ((stress(1)-avg)^2+stress(3)^2)^(1/2);
     %p1  = avg+R;
     %p2 = avg-R;
     
     
     vonM = sqrt(stress(1)^2  +   stress(2)^2 -  stress(1)*stress(2)  +   3*(stress(3))^2);
     % vonM2 = sqrt(p1^2-p1*p2+p2^2);
     
     % note: vonM and vonM2 should be the same. 
     vonM_stored(e) = vonM;    
     
    % ---------------------------------------
    % plot the element outline and the displacments
    % ---------------------------------------
     hold on
     coordD = zeros(5,1);   
     for temp = 1:4
        coordD(temp,1) =  coord(temp,1) + multiplierScale*U(2*arrayCoordNumber(temp)-1); % X value
        coordD(temp,2) =  coord(temp,2) + multiplierScale*U(2*arrayCoordNumber(temp)); % Y value
     end    
     
     coord2 = coord;
     coordD(5,:) = coordD(1,:) ;
     coord2(5,:) = coord2(1,:); 
     plot(coord2(:,1),coord2(:,2),'-g');  
     plot(coordD(:,1),coordD(:,2), '-b');   
     
end
axis equal
tti= strcat('Displacement of the elements shown in Blue');
title(tti);
hold off 

if(1==0)

    x = elemcenterLocations(:,1)';
    y = elemcenterLocations(:,2)';
    z = vonM_stored'; 

    figure(2)
    scatter(x,y,[],z)
    xlim([0,L])
    ylim([0,h])
    axis equal

    figure(3)
    % http://www.mathworks.com/matlabcentral/fileexchange/38858-contour-plot-for-scattered-data
    tri=delaunay(x,y);           % triangulate scattered data
    vonMax = max(vonM_stored);
    step = floor(vonMax/30);
     v=0:step:vonMax; % contour levels
    [C,h]=tricontour(tri,x,y,z,v);
    %clabel(C,h)
    xlim([0,L])
    %ylim([0,h])
    axis equal

    % ------------------------------------------------------------------------
    % -----------------Final calculations ------------------------------------
    % ------------------------------------------------------------------------
    % ----------------------------------------
    disp('Max Von Mises Stress is')
    max(vonM_stored)

    if(problem ==1)  

        % Displacement at L
        disp('Displacement at L')
        U(2*nely), U(2*(nelx+1))

        % ----------------------------------------
        % Displacement at d
        % ----------------------------------------
        %
        % Solve for j
        % d = (j-1)*L/(nelx);
        %
        % d*nelx/L+1 = j

        j =  d*nelx/L+1 ;
        left = floor(j);
        right = ceil(j);

        averageXDisplacement = (U(left*2-1)+U(right*2-1))/2;
        averageYDisplacement = (U(left*2)+U(right*2))/2;

        disp('Displacement at d is , (x,y)');
        [averageXDisplacement, averageYDisplacement] %#ok<NOPTS>

        % ----------------------------------------
        % Displacement at a
        % ----------------------------------------

        % Solve for a instead of d
        j =  a*nelx/L+1 ;
        left = floor(j);
        right = ceil(j);

        averageXDisplacement = (U(left*2-1)+U(right*2-1))/2;
        averageYDisplacement = (U(left*2)+U(right*2))/2;

        disp('Displacement at a is , (x,y)');
        [averageXDisplacement averageYDisplacement] %#ok<NOPTS>

    else

        % ------------------------------------------------------
        % Problem 2
        %
        % The A node calculation doesn't work for the mesh refinement, it is
        % hard coded becasue I was having problems. 
        % ------------------------------------------------------
        Lnode = 0;
        Anode = 9;
        Dnode = 0;
        for i = 1:nn
            % get the xy locations as a row
            xy = globalPosition(i,:);
            xLoc = xy(1); yLoc = xy(2);

            % Find L
            if( yLoc  == 0 && xLoc ==L )
                Lnode = i;
            end

            % Find d
         if( yLoc  == 0 && xLoc ==d )
            Dnode = i;
         end

             % Find a, inside the notch
            if( yLoc  == 0.1670 && xLoc ==4 )
                Anode = i;
            end
        end

        % Displacement at L
        disp('Displacement at L')
        [U(2*Lnode-1) U(2*(Lnode))]  %#ok<NOPTS>

        % ----------------------------------------
        % Displacement at d
        % ----------------------------------------

        disp('Displacement at d is , (x,y)');
         [U(2*Dnode-1) U(2*Dnode)]  %#ok<NOPTS>

        % ----------------------------------------
        % Displacement at a
        % ----------------------------------------

          disp('Displacement at a is , (x,y)');
        [U(2*Anode-1) U(2*Anode)]  %#ok<NOPTS>
    end
end