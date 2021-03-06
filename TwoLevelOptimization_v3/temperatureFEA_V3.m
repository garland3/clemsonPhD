function [T]  = temperatureFEA_V3(DV, config, matProp,loop, loadcase)


u0 =0; % value at essentail boundaries
nn = (config.nelx+1)*(config.nely+1); % number of nodes
ne = config.nelx*config.nely; % number of elements

% ------------------------------------
% -- Make a mapping between elements and global matrixes
% ------------------------------------
% IEN holds the node numbers for each element. 
% Each row is a new element
% The first column is element 1's global node number
% Second column is elemnt 2's global node number
%  and ....

% count = 1;
% elementsInRow = settings.nelx+1;
% IEN = zeros(nn,4);
% % Each row, so nely # of row
% for i = 1:settings.nely
%      rowMultiplier = i-1;
%     % Each column, so nelx # of row
%     for j= 1:settings.nelx        
%         IEN(count,:)=[rowMultiplier*elementsInRow+j,
%                       rowMultiplier*elementsInRow+j+1,
%                       (rowMultiplier +1)*elementsInRow+j+1,
%                        (rowMultiplier +1)*elementsInRow+j];
%         count = count+1;
%     end
% end
% 
% numNodesInRow = settings.nelx+1;
% numNodesInColumn = settings.nely+1;
% XLocations=zeros(numNodesInRow,numNodesInColumn);
% YLocations=zeros(numNodesInRow,numNodesInColumn);

% Find and store the global positions of each node
% Each element square is 1 by 1 units
% Store both the X and Y positions
% globalPosition = zeros(nn,2);
% count = 1;
% for i = 1:numNodesInColumn  % y
%     for j= 1:numNodesInRow % x
%         globalPosition(count,:) = [j-1 i-1];
%         count = count+1;        
%         XLocations(j,i) = j-1;
%         YLocations(j,i) = i-1;
%     end
% end 

% Specifiy the constrained nodes where there are essential boundary
% conditions

K = zeros(nn,nn);
row = config.nelx+1;
column = config.nely+1;

  loadingScenario = 'sourceEverywhereSinkMiddle';
%  loadingScenario = 'pressure';
% 'topAndRightLoadsSinkBottomLeft'
% 'sourceMiddleRightSinkMiddleLeft'

if strcmp(loadingScenario,'heatMiddleSinksCorners')
    % set a heat source in the middle and sinks on the four corners. 
    F = zeros(nn,1);
    F([ceil(row/2)+(ceil(column/2)*row) (ceil(row/2)+1)+(ceil(column/2)*row)]) =  20; % heat source in the middle
   
    % Set the 4 corners. 
    quartY = ceil(config.nely/4);
    Essential =   [1 row (column-1)*row+1 column*row] ; % the 4 corners
    
    
elseif (strcmp(loadingScenario,'pressure'))
  
    
    middleX = config.nelx/2;
    middleY = config.nely/2;
    radius = 10;
    error=1;
    doPlotHeatPositions=0;
    EssentialCold=[];
    EssentialHot=[];
    
    if(doPlotHeatPositions==1)
        heatPositions = zeros(size(DV.x));
        sinkPositions = zeros(config.nely+1,config.nelx+1);
    end
    
    ForceValue=0.01;
      F = zeros(nn,1);
    for ii = 1:config.nelx
        for jj = 1:config.nely
            distanceFromMiddle=sqrt((ii-middleX)^2+(jj-middleY)^2);
            if(distanceFromMiddle<radius+error)
                xNodeNum = (jj-1)*row+ii-1;
                 F(xNodeNum)=ForceValue;
% EssentialHot=[EssentialHot xNodeNum];
                
                heatPositions(jj,ii)=ForceValue;
            end
        end
    end
    
    tt=   1:row :row*(column); % ... % left side
    t2=tt-1;
    EssentialCold=[tt t2 ];
EssentialCold=[ EssentialCold 1: row (column-1)*row: nn];

    EssentialCold=EssentialCold(EssentialCold>0);
      EssentialCold=EssentialCold(EssentialCold<=nn);
    EssentialCold = unique(EssentialCold);


Essential=[ EssentialCold EssentialHot];
%     Essential=Essential(Essential>0);
%       Essential=Essential(Essential<=nn);
    Essential = unique(Essential);




    
    if(doPlotHeatPositions==1)
        subplot(1,2,1);
        imagesc(heatPositions);
        
        for ii = 1:config.nelx+1
            for jj = 1:config.nely+1
                indexNumber = (jj-1)*row+ii;
                if(  any(indexNumber==Essential)==1)
                    sinkPositions(jj,ii)=1;
                end
            end
        end
      
          subplot(1,2,2);
         imagesc(sinkPositions);
    end
    
    

elseif (strcmp(loadingScenario,'topAndRightLoadsSinkBottomLeft'))
    
    
   % ----------------------
   % Heat source in top middle and right middle. Sink in bototm left
   % ----------------------
      F = zeros(nn,1)*0.02;
      F([row+(ceil(column/2)*row) ]) =  0.5; % heat source in the middle right
      F([(column-1)*row+ceil(row/2) ]) =  0.5; % heat source in the top middle

    Essential = [1 2 row+1];
    
elseif  (strcmp(loadingScenario,'sourceMiddleRightSinkMiddleLeft'))
    
       % ----------------------
       % Heat source in right middle. Sink in left middle
       % ----------------------
      F = zeros(nn,1)*0.02;
       F([row+(ceil(column/2)*row) ]) =  1; % heat source in the middle right
      

     Essential = [1+(ceil(column/2)*row)]; % sink in left middle
elseif (strcmp(loadingScenario,'sourceMiddleSinkBoundaries'))
     F = zeros(nn,1);
    F([ceil(row/2)+(ceil(column/2)*row) (ceil(row/2)+1)+(ceil(column/2)*row)]) =  20; % heat source in the middle   
    % Set the 4 corners. 
    bottomRow = 1:row;
    topRow = (column-1)*row+1:column*row;
    left = 1:row:(column-1)*row+1;
    right = row:row:column*row;
    Essential =   [bottomRow topRow left right] ; % 4 sides    
elseif (strcmp(loadingScenario,'sourceEverywhereSinkMiddle'))    
    F = ones(nn,1)*0.001;
    Essential=[ceil(row/2)+(ceil(column/2)*row) (ceil(row/2)+1)+(ceil(column/2)*row)]; % heat sink in the middle   
    
elseif (strcmp(loadingScenario,'sourceEverywhereSinkBottomLeft'))  
%      F = ones(nn,1)*0.001;
        F = ones(nn,1)*0.1;
      Essential=[1 2 row+1];
end

Essential = unique(Essential);
alldofs     = [1:nn];
Free    = setdiff(alldofs,Essential);
%F = ones(nn,1); % add a constant heat source everywhere


  
% xLoc = 1;
% yLoc = 1;
% loop over the elements
for e = 1:ne
    
      % Get the precalculated element stiffness matrix. 
     % ke = elementK_heat();
      [elx,ely]= DV.GivenNodeNumberGetXY(e);
      
        %xx= nelx; yy = nely;
      if(config.doUseMultiElePerDV==1) % if elements per design var. 
         [elx,ely] = DV.GetDesignVarPositionGivenXYElement(config,elx,ely);
      end
      
      ke = matProp.effectiveHeatKEmatrix(  DV.w(ely,elx), config);

      
      
      % Insert the element stiffness matrix into the global.    
      node = DV.IEN(e,:);
      
      % for the x location
      % The first number is the row - "y value"
      % The second number is the column "x value"
       K(node,node) = K(node,node) + DV.x(ely,elx)^config.penal*ke;
       
%        xLoc = xLoc+1;
%        if(xLoc>settings.nelx)
%            xLoc = 1;
%            yLoc = yLoc+1;
%        end

end
 
K = sparse(K);
F = sparse(F);

K_ff = K(Free,Free);
% K_fe = K(Free,Essential);
F_f = F(Free);
  
% if(settings.useGPU ==1)
%     % GPU matrix solve. 
%     K_ff_gpu = gpuArray(K_ff);
%     F_f_gpu = gpuArray(F_f);
%     T_gpu = K_ff_gpu\F_f_gpu;
%     T(Free) = gather(T_gpu);
% else
    % normal matrix solve
     T(Free) = K_ff \ F_f;

% end


% if (strcmp(loadingScenario,'pressure'))
%     
%     T(EssentialCold) = u0;
%     T(EssentialHot) = 100;
%     
% else
    T(Essential) = u0;
% end


  
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

% TcontourMatrix = zeros(numNodesInRow,numNodesInColumn);
% 
%  
%   
% for j = 1:numNodesInColumn % y
%       rowMultiplier = j-1;
%      for i = 1:numNodesInRow % x
%          nodeNumber = i+numNodesInRow*rowMultiplier;
%          TcontourMatrix(i,j) = T(nodeNumber);
%      
%      end
% end

% if(settings.doPlotHeat==1 && mod(loop,settings.iterationsPerPlot) ==0)
%      subplot(1,2,1)
%      % plot the coutour graph
%      contour(XLocations,YLocations,TcontourMatrix);
%      tti= strcat('Heat contours. Number of elements=', int2str(ne));
%      title(tti);
%      
%     % 
%     % plot the surf graph
% %     subplot(1,2,2)
% %     surf(XLocations,YLocations,TcontourMatrix);
% end
% subplot(1,2,1);
%  surf(XLocations,YLocations,TcontourMatrix);
% tti= strcat('Heat surface.  Number of elements =', int2str(ne));
% title(tti);
% 
% diary off
T = transpose(T);
% T = transpose(T);
