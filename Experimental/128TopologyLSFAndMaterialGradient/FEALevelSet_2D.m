% Copyright Anthony Garland 2015
% ------------------------
function [U, g1_local_square,g2_local_square, volFracV1, volFracV2] = FEALevelSet_2D(structure,lsf, volFracArray,  doplot, ...
    alphaPenalty,beta, countMainLoop, mode, resolutionMultiplier)
% structure - shows the structure's boundary by a binary relationship, 0 = void, 1 = material
% lsf - is the level set function. I need this level set function so that I
% can calcualte the normal direction and mean curvature
% volFracArray - tells the volume fraction at each element in the domain
% doplot - indicates if we should plot or not
% alphaPenalty - indicates the strength of the alpha penalty function which
% penalizes the volume fraction to encourage smoothness

% --------------------------
% Plotting information
% --------------------------

% 2D grid with a 3D function 
% http://www.mathworks.com/help/matlab/ref/interp2.html
recvid = 1; % Record a video of the figure 1, record view, 1 = yes

iterationsPerPlot = 10;

doplotDisplacement = doplot; % Set to 1 to show. Runs much slower
plotStress = doplot; % Set to 1 to plot the stress graphs
%plotStress = 1; % Set to 1 to plot the stress graphs
plotElasticMod = doplot;
%plotElasticMod = 1;
plotVolFraction = doplot;
%plotVolFraction = 1; % override
plotStructure = doplot;
%plotStructure = 1; % override
plotLSF = doplot; 
%plotLSF = 1; % override
plotNormalDirection = doplot;

plotSensitivity = doplot;

% if recvid==1
%     vidObj = VideoWriter('results.avi');    %Prepare the new file for video
%     vidObj.FrameRate = 5;
%     vidObj.Quality = 100;
%     open(vidObj);
%     vid=1;
% end

rM = resolutionMultiplier;
subplotY = 2; % Suplot matrix setup
subplotX = 2;
subplotCount = 1;
multiplierScale = 10; % When plotting the diplaced elements, exaggerate the displacements by this scale

if (mode ==1) % optimize material only
    plotStress = 1; % Set to 1 to plot the stress graphs
    plotElasticMod = 1;
    plotVolFraction = 1; % override
elseif( mode ==2) % optimize topology only
    plotStructure = 1;
    plotLSF = 1;
    plotStress = 1; % Set to 1 to plot the stress graphs
    plotSensitivity = 1; % override
elseif (mode ==3)
    plotStructure = 1;
    plotLSF = 1;
    plotStress = 1;
     plotVolFraction = 1; % override
    
    
    
end

[nely,nelx] = size(volFracArray);
% nely = nely/rM;
% nelx = nelx/rM;

h = nely;
L = nelx;

if(doplot ==1)
    figure(1)
    clf    
end
imageXaxis = [0 nelx];
imageYaxis = [0 nely];

% ------------------------------
% FEA, Beam material properties and physical properties
% use SI units
% ------------------------------
t = 5; % mm, thickness of the beam
FappliedLoad = 200; % N. 
Enylon = 1700; % elastic mod of nylon, Mpa
Epla =  3368; % elastic mod of pla, MPa
E_empty = 1; % elastic mod for the region with no material to prevent sigularity of the FEA
v = 0.3; % Assume 0.3 for Poisson's ratio for both materials

% ------------------------
% Calculate the derivative of the D matrix with respect to a material
% change. I call this D, but in the paper, they call it A. 
% ------------------------
D_dd =  [ 1 v 0;
              v 1 0;
              0 0 1/2*(1-v)]*(Epla-Enylon)/(1-v^2);
          
Aomega = D_dd;


% Calculate the Laplace of the vol Fraction composition. (needed for the G1
% term calculation)
stepSize = 1;
LaplaceVolFract = del2(volFracArray,stepSize);

% Calculate the gradient as well ofr the volume fraction, this is neededd
% for the G2 term)
[gvolFracX, gvolFracY] = gradient(volFracArray,stepSize);
volFracGradientSquared = gvolFracX.^2+gvolFracY.^2;


% Find the normal direction of the lsf,  = nabla(lsf)/abs(nabla(lsf))
[gx_lsf, gy_lsf] = gradient(lsf);
denominator_g_lsf = (gx_lsf.^2 +gy_lsf.^2).^(1/2);
gx_lsf = gx_lsf./denominator_g_lsf;
gy_lsf = gy_lsf./denominator_g_lsf;

curvature_lsf = divergence(gx_lsf,gy_lsf); % calculate the divergence of the normal direction to get the curvature
curvature_lsf(isnan(curvature_lsf)) = 0 ; % remove the NaN
curvature_lsf = curvature_lsf(2:end-1,2:end-1);

curvature_lsf = curvature_lsf*0+1;

% if (plotNormalDirection == 1)
%      figure(1)
%     subplot(subplotY,subplotX,subplotCount); subplotCount=subplotCount+1; 
%     quiver(gx_lsf, gy_lsf);
%     title('Normal direction of lsf');
% end
    

% Handle the mapping between nodes and elements and node
% locations

% nelx; %  Number of elements in the x direction
% nely; %  Number of elements in the y direction
nn = (nelx+1)*(nely+1); % number of nodes
ne = nelx*nely; % number of elements 

% ---------------
% Calcualte the Elastic mod at each element
% Find the volume fractions of each material
% ---------------
E_atElement = volFracArray;
volFracV1= 0;
volFracV2 = 0;

 for i = 1:nelx
    for j = 1:nely           
         structureLocal = structure(j,i);
            if(structureLocal == 0) % if void region
                E_atElement(j,i) = E_empty;                     
            else % if a filled region
                  volFraclocal = volFracArray(j,i);
                  volFracV1 = volFracV1 +volFraclocal; % sum up the total use of material 1 (PLA)
                  volFracV2 = volFracV2 + (1- volFraclocal); % sum up the total use of material 2 (Nylon)
                 E_atElement(j,i)= Epla*volFraclocal+(1-volFraclocal)*Enylon;  % simple mixture ratio 
            end     
    end
 end    
 
 % normalize the volume fraction by the number of elements
 volFracV1 = volFracV1/ne;
 volFracV2 = volFracV2/ne;
 
 if(plotLSF ==1 && mod(countMainLoop,iterationsPerPlot) ==0 )
    figure(1)
    subplot(subplotX,subplotY,subplotCount); subplotCount=subplotCount+1;
   imagesc(lsf); axis equal; axis tight; axis off;
   colormap winter
    colorbar 
   % caxis([-1 1 ]);
   title('LSF')
end

% --------------------------
% Plot Elastic mod
% --------------------------
if(plotElasticMod ==1 && mod(countMainLoop,iterationsPerPlot) ==0)   
    figure(1)
    subplot(subplotY,subplotX,subplotCount); subplotCount=subplotCount+1;       
    imagesc(imageXaxis,imageYaxis,E_atElement);      

        % Set the X and Y axis labels
        xlabel('X','FontSize',10, 'FontName','Arial')
        ylabel('Y','FontSize',10, 'FontName','Arial')
         
         % http://www.mathworks.com/help/matlab/ref/colormap.html    
         colormap winter
%           cmap = colormap;
%           cmap(1,:) = 1;
%           colormap(cmap);
          % freezeColors
          caxis([Enylon Epla ])
         set(gca,'YDir','normal'); % http://www.mathworks.com/matlabcentral/answers/94170-how-can-i-reverse-the-y-axis-when-i-use-the-image-or-imagesc-function-to-display-an-image-in-matlab
          axis equal 
         colorbar
         title('Elastic mod over the domain')
          xlim([0,L])
         ylim([0,h])
end

% --------------------------
% Plot structure
% --------------------------

if(plotStructure ==1 && mod(countMainLoop,iterationsPerPlot) ==0)
    figure(1)
    subplot(subplotX,subplotY,subplotCount); subplotCount=subplotCount+1;
    imagesc(structure); axis equal; axis tight; axis off;
      set(gca,'YDir','normal');
       colorbar;
    title('Structure')
end

% --------------------------
% Plot Volume Fraction
% --------------------------
if(plotVolFraction ==1 && mod(countMainLoop,iterationsPerPlot) ==0)  
    figure(1)
    subplot(subplotY,subplotX,subplotCount); subplotCount=subplotCount+1;       
    imagesc(imageXaxis,imageYaxis,volFracArray);
    set(gca,'YDir','normal'); % http://www.mathworks.com/matlabcentral/answers/94170-how-can-i-reverse-the-y-axis-when-i-use-the-image-or-imagesc-function-to-display-an-image-in-matlab

    % Set the values in in ches and give mm label   
    ax1 = gca;
    set(ax1, 'XLim', [0 nelx ],'YLim', [0 nely] );  
    % Set the X and Y axis labels
    xlabel('X (mm)','FontSize',10, 'FontName','Arial')
    ylabel('Y (mm)','FontSize',10, 'FontName','Arial')
    title('Volume Fraction Composition')
     caxis([0 1 ])
     colorbar
    xlim([0,L])
    ylim([0,h])
    axis equal 
    
end

% ------------------------------------
% -- Make a mapping between elements and global matrixes
% ------------------------------------

% ElemToNodeMap holds the node numbers for each element. 
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

E_atElementsArray = zeros(ne,1);
LaplaceVolFrac_atElementsArray = zeros(ne,1);
GradientSquredVolFact_atElementsArray = zeros(ne,1);
curvature_lsf_atElements = zeros(ne,1);


count = 1;
numNodesInRow = nelx+1;
numNodesInColumn = nely+1;
ElemToNodeMap = zeros(nn,4); % Element To Node Mapping
% Each row, so nely # of row
for i = 1:nely
     rowMultiplier = i-1;
    % Each column, so nelx # of row
    for j= 1:nelx              
        % Store the E value for this element
        E_atElementsArray(count) = E_atElement(i,j);    
        LaplaceVolFrac_atElementsArray(count) = LaplaceVolFract(i,j); % store the laplace of the vol fraction in a single column matrix
        GradientSquredVolFact_atElementsArray(count) = volFracGradientSquared(i,j);
        curvature_lsf_atElements(count) = curvature_lsf(i,j);
         % Store node mapping
        ElemToNodeMap(count,:)=[rowMultiplier*numNodesInRow+j, ...
                      rowMultiplier*numNodesInRow+j+1, ...
                      (rowMultiplier +1)*numNodesInRow+j+1, ...
                       (rowMultiplier +1)*numNodesInRow+j];
        count = count+1;
    end
end

% Find and store the global positions of each node
% Each element rectangle. The aspect ratio square for each element
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
   
% ---------------------------------------------
% Initialize the global force, k, and B_stored matrixes
% ---------------------------------------------
ndof = nn*2; % Number of degrees of freedome. 2 per node. 
F2 = zeros(ndof,1);
K = zeros(ndof,ndof);
 
% ---------------------------------------------
% Find the essential boundary conditions
% Bridge problem
% ---------------------------------------------
u0 = 0; % value at essential boundaries
row = nelx+1;
Essential2 = [1 2]; % bottom left corner is fixed
Essential2 = [Essential2 row*2 row*2-1]; % bottom right corner is only fixed in the y direction
F2( (floor(row/2)+1)*2) = -FappliedLoad; % force down in the bottom middle

% for i = 1:nn
%     % get the xy locations as a row
%     xy = globalPosition(i,:);
%     xLoc = xy(1); yLoc = xy(2);
%     
%     if( xLoc ==0)
%         % exists returns 1 if Essential2 is a varriable
%         temp = exist('Essential2', 'var');
%         
%         if(temp==1)
%             % only constrain in the X direction, 2*i-1
%              Essential2 = [Essential2  2*i-1]; %#ok<AGROW>
%         else
%             Essential2 =[2*i-1];
%         end
%     end
%     
%     % if in the left top, then add the downward force
%     if(yLoc == 2 && xLoc == 0)
%          F2(i*2) = -FappliedLoad;
%     end
%     
%     % if bottom right, the add a boundary condition to prevent up-down
%     % movement
%     if(yLoc == 0 && xLoc == 6)
%          Essential2 = [Essential2 (i*2)];
%     end
%         
% end
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
  
% loop over the elements
for e = 1:ne         
      % loop over local node numbers to get their node global node numbers
      coord = zeros(4,2);
      for j = 1:4
          % Get the node number
          nodeNumber = ElemToNodeMap(e,j);
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
         
      % find the E for this element. This was precalculated earlier
      E = E_atElementsArray(e);
      D = [ 1 v 0;
              v 1 0;
              0 0 1/2*(1-v)]*E/(1-v^2);
          
      % Loop over the guass points
      ke = zeros(8,8);
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
        
%        % Calculate the Area as determinate(Jacobian)*4
%        area =  det(J)*4;
%        volume = area*t; % t is the thickness
%        bodyForce = volume*density;
%        f_element = zeros(8,1);
%        f_element(2:2:8,1) = -bodyForce/4;
      
     
     % Insert the element stiffness matrix into the global.    
     nodes1 = ElemToNodeMap(e,:);
     xDOF = nodes1*2-1;
     yDOF = nodes1*2;
    
     % I cannot use the union, or else the order get messed up. The order
     % is important. Same in the actual topology code when you are
     % calculating the objectiv
      dofNumbers = [xDOF(1) yDOF(1) xDOF(2) yDOF(2) xDOF(3) yDOF(3) xDOF(4) yDOF(4)];
      
      % multiply t*ke and add to global matrix. t = thickness or t_z
      K(dofNumbers,dofNumbers) = K(dofNumbers,dofNumbers) + t*ke;
      
      % Don't hadd the body force right now
      % F2(dofNumbers,1) = F2(dofNumbers,1) +f_element;        
end
  
K_ff = K(Free,Free);
% K_fe = K(Free,Essential2);
F_f = F2(Free);
  
U(Free) = K_ff \ F_f; % solve the FEA
U(Essential2) = u0;  

% Max displacement can be used a preformance measure in some cases. 
% maxDisplacement1 = max(max(U));
% minDisplacement1 = min(min(U));
% maxDisplacement = max([maxDisplacement1 -minDisplacement1]); % get the max in the positive or negative direction of U displacements

stress_stored = zeros(ne,3);
vonM_stored = zeros(ne,1);
g1_localstored = zeros(ne,1);
g2_localstored = zeros(ne,1);
elemcenterLocations = zeros(ne,2);

if(doplotDisplacement ==1 && mod(countMainLoop,iterationsPerPlot) ==0)
    figure(1)
    subplot(subplotY,subplotX,subplotCount); subplotCount=subplotCount+1;
end

strainEnergy = 0;

% --------------------------------
% Post processing
% Calculate stress, vonMises, and sensitivies 
% --------------------------------

% loop over the elements
 for e = 1:ne
    arrayCoordNumber = zeros(1,4);
     xsum = 0;
     ysum = 0;
     % loop over local node numbers to get their node global node numbers
     local_u = zeros(4,1);
     for j = 1:4
         % Get the node number
         nodeNumber = ElemToNodeMap(e,j);
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
     nodes1 = ElemToNodeMap(e,:);
     xDOF = nodes1*2-1;
     yDOF = nodes1*2;
     dofNumbers = [xDOF(1) yDOF(1) xDOF(2) yDOF(2) xDOF(3) yDOF(3) xDOF(4) yDOF(4)];
     
     d_local = U(dofNumbers);
     
     % Compute D
     E = E_atElementsArray(e);
     D = [ 1 v 0;
              v 1 0;
              0 0 1/2*(1-v)]*E/(1-v^2);
     
     % Strain is B*d, NOTE: d is transposed
     strain = B*d_local';         
     stress = D*strain;
     vonM = sqrt(stress(1)^2  +   stress(2)^2 -  stress(1)*stress(2)  +   3*(stress(3))^2);
     vonM_stored(e) = vonM;    
     
     % G1_local
     g1_local = alphaPenalty*LaplaceVolFrac_atElementsArray(e)+ strain'*(Aomega*(strain));
     g1_localstored(e) = g1_local;     
    
     
     localStrainE = strain'*stress;
     strainEnergy = strainEnergy+localStrainE;
     
      % G2_local
     g2_local = localStrainE-alphaPenalty*GradientSquredVolFact_atElementsArray(e) - beta*curvature_lsf_atElements(e);
     g2_localstored(e) = g2_local;
             
    % Store the transpose, to make things work nice. 
    stress_stored(e,:) = stress';
    
    % ---------------------------------------
    % plot the element outline and the displacments
    % ---------------------------------------    
    if(doplotDisplacement ==1 && mod(countMainLoop,iterationsPerPlot) ==0)
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
     
 end
 
if(doplotDisplacement ==1 && mod(countMainLoop,iterationsPerPlot) ==0)
    axis equal
    tti= strcat('Displacement of the elements shown in Blue');
    title(tti);
    hold off 
end





% Show final results
if(plotStress ==1 && mod(countMainLoop,iterationsPerPlot) ==0)

     x = elemcenterLocations(:,1)';
     y = elemcenterLocations(:,2)';
     z = vonM_stored'; 
% 
%     subplot(subplotY,subplotX,subplotCount); subplotCount=subplotCount+1;
%     scatter(x,y,[],z)
%     xlim([0,L])
%     ylim([0,h])
%     axis equal

    subplot(subplotY,subplotX,subplotCount); subplotCount=subplotCount+1;
  %  http://www.mathworks.com/matlabcentral/fileexchange/38858-contour-plot-for-scattered-data
        tri=delaunay(x,y);           % triangulate scattered data
        vonMax = max(vonM_stored);
        stepX = 2;
         v=0:stepX:vonMax; % contour levels
        [C,hh]=tricontour(tri,x,y,z,v);
        %clabel(C,h)
        xlim([0,L])
        ylim([0,h])
        axis equal

   
end

% -----------------------
% Convert the g1_local into a matrix rather than just a long column
% -----------------------
g1_local_square = zeros(nely, nelx);
g2_local_square = zeros(nely,nelx);
count = 1;
for i = 1:nely   
    for j= 1:nelx  
        g1_local_square(i,j) = g1_localstored(count);
        g2_local_square(i,j) = g2_localstored(count);
        count = count+1;
    end
end


if(plotSensitivity ==1 && mod(countMainLoop,iterationsPerPlot) ==0)   
    figure(1)
    subplot(subplotY,subplotX,subplotCount); subplotCount=subplotCount+1;       
    imagesc(imageXaxis,imageYaxis,g2_local_square);      

        % Set the X and Y axis labels
        xlabel('X','FontSize',10, 'FontName','Arial')
        ylabel('Y','FontSize',10, 'FontName','Arial')
         
         % http://www.mathworks.com/help/matlab/ref/colormap.html    
         colormap winter
          %cmap = colormap;
          %cmap(1,:) = 1;
          %colormap(cmap);
          % freezeColors
          caxis([min(min(g2_local_square)) max(max(g2_local_square)) ])
         set(gca,'YDir','normal'); % http://www.mathworks.com/matlabcentral/answers/94170-how-can-i-reverse-the-y-axis-when-i-use-the-image-or-imagesc-function-to-display-an-image-in-matlab
          axis equal 
         colorbar
         title('Sensitivity')
          xlim([0,L])
         ylim([0,h])
end

%  if ( recvid==1 &&  mod(countMainLoop,iterationsPerPlot) ==0)
%         F(vid) = getframe(figure(1)); %#ok<AGROW> %Get frame of the topology in each iteration
%         writeVideo(vidObj,F(vid)); %Save the topology in the video
%         vid=vid+1;
%  end
%  
% if(countMainLoop >=499 &&  recvid==1) 
% 
%     close(vidObj);  %close video
% end


end
