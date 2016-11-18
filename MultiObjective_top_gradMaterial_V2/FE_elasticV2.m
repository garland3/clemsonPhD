function [T, maxF,maxT]=FE_elasticV2(designVars, settings, matProp, loadingCase)

u0 =0; % value at essentail boundaries
nn = (settings.nelx+1)*(settings.nely+1); % number of nodes
ne = settings.nelx*settings.nely; % number of elements
ndof = nn*2; % Number of degrees of freedome. 2 per node.

% Specifiy the constrained nodes where there are essential boundary
% conditions
F = zeros(ndof,1);
K = zeros(ndof,ndof);
row = settings.nelx+1;
column= settings.nely+1;

bottomFixed=0;
% loading condition
if loadingCase == 111
    FappliedLoad = 500;
    tt=   1:2*row :2*row*(column); % ... % left side
    t2=tt+1;      
    t3 = [];
    t4 = [];
    u3 = 0;
    Essential=[tt t2 t3 t4];
    Essential = unique(Essential);
    % Down in the top right corner
    F(ndof) = FappliedLoad;
    
elseif loadingCase ==112
    FappliedLoad = -500;
    tt=   1:2*row :2*row*(column); % ... % left side
    t2=tt+1;
    Essential=[tt t2];
    Essential = unique(Essential);    
    % down on bottom right
    F(row*2) = FappliedLoad;
    
elseif loadingCase ==113
    FappliedLoad = 500;
    tt=   1:2*row :2*row*(column); % ... % left side
    t2=tt+1;
    Essential=[tt t2];
    Essential = unique(Essential);    
    F( (row*2)*floor(column/2)) = FappliedLoad; % middle to the right  
    

    % -------------------------------------------------
    %
    % 300s are the shoe loading case. 
    %
    % -------------------------------------------------
elseif loadingCase ==300
      FappliedLoad = -500;   
      bottomFixed = 1;
     
      fnodes =  (row*(column-1)+1)*2:2:(row*(column-1)+floor(row/4))*2; % top left quarter, y degrees of freedom, pointing down. 
      count = size(fnodes,2);
       F(fnodes) =FappliedLoad/count;
      
elseif loadingCase ==301
        FappliedLoad = -500;   
      bottomFixed = 1;
     
      fnodes =   (row*(column-1)+1)*2:2:(row*(column-1)+floor(row/2))*2; % top left half, y degrees of freedom, pointing down. 
      count = size(fnodes,2);
       F(fnodes) =FappliedLoad/count;
elseif loadingCase ==302
      FappliedLoad = -500;   
      bottomFixed = 1;
     
     fnodes= (row*(column-1)+1)*2:2:ndof; % whole top, y degrees of freedom, pointing down. 
      count = size(fnodes,2);
       F(fnodes) =FappliedLoad/count;
elseif loadingCase ==303
         FappliedLoad = -500;   
      bottomFixed = 1;
     
      fnodes =  (row*(column-1)+1+floor(row/2))*2:2:ndof; % top right half, y degrees of freedom, pointing down. 
      count = size(fnodes,2);
       F(fnodes) =FappliedLoad/count;
elseif loadingCase ==304
     FappliedLoad = -500;   
      bottomFixed = 1;
     
      fnodes =  (row*(column-1)+1+floor(3*row/4))*2:2:ndof; % top right quarter, x and y degrees of freedom, pointing down and left. 
      count = size(fnodes,2);
       F(fnodes) =FappliedLoad/count;
       
elseif loadingCase ==305
     FappliedLoad = -500;   
      bottomFixed = 1;
     
      fnodes =  (row*(column-1)+floor(row/2))*2:2:ndof; % top right half, x  degree of freedome only, pushing off. 
      count = size(fnodes,2);
       F(fnodes) =FappliedLoad/count;
       
    % -------------------------------------------------
    %
    % 400s are the BRIDGE
    %
    % -------------------------------------------------
elseif loadingCase ==400
    
   
     FappliedLoad = -500;   
      bottomFixed = 1;
      
        start = (row*(floor(column/2))+1)*2;
        endnode=(row*(1+floor(column/2)))*2;
        
        startOffset = 0;
        endOfset = -(floor(3*row/4))*2;
        
     fnodes= start+startOffset :2:endnode+endOfset ; % middle, 1/4 load down
     count = size(fnodes,2);
     F(fnodes) =FappliedLoad/count;

elseif loadingCase ==401
     
     FappliedLoad = -500;   
      bottomFixed = 1;
      
        start = (row*(floor(column/2))+1)*2;
        endnode=(row*(1+floor(column/2)))*2;
        
        startOffset =(floor(1*row/4))*2;
        endOfset = -(floor(2*row/4))*2;
        
     fnodes= start+startOffset :2:endnode+endOfset ; % middle, 1/4 load down
     count = size(fnodes,2);
     F(fnodes) =FappliedLoad/count;
elseif loadingCase ==402
      FappliedLoad = -500;   
      bottomFixed = 1;
      
        start = (row*(floor(column/2))+1)*2;
        endnode=(row*(1+floor(column/2)))*2;
        
        startOffset =(floor(2*row/4))*2;
        endOfset = -(floor(1*row/4))*2;
        
     fnodes= start+startOffset :2:endnode+endOfset ; % middle, 1/4 load down
     count = size(fnodes,2);
     F(fnodes) =FappliedLoad/count;
      
elseif loadingCase ==403
      FappliedLoad = -500;   
      bottomFixed = 1;
      
        start = (row*(floor(column/2))+1)*2;
        endnode=(row*(1+floor(column/2)))*2;
        
        startOffset =(floor(3*row/4))*2;
        endOfset = 0;
        
     fnodes= start+startOffset :2:endnode+endOfset ; % middle, 1/4 load down
     count = size(fnodes,2);
     F(fnodes) =FappliedLoad/count;
     
 elseif loadingCase ==404
      FappliedLoad = -500;   
      bottomFixed = 1;
      
        start = (row*(floor(column/2))+1)*2;
        endnode=(row*(1+floor(column/2)))*2;
        
        startOffset =-11;
        endOfset = 0;
        
     fnodes= start+startOffset :2:endnode+endOfset ; % middle, 1/4 load down
     count = size(fnodes,2);
     F(fnodes) =FappliedLoad/count;
  elseif loadingCase ==405
      FappliedLoad = 500;   
      bottomFixed = 1;
      
        start = (row*(floor(column/2))+1)*2;
        endnode=(row*(1+floor(column/2)))*2;
        
        startOffset =-11;
        endOfset = 0;
        
     fnodes= start+startOffset :2:endnode+endOfset ; % middle, 1/4 load down
     count = size(fnodes,2);
     F(fnodes) =FappliedLoad/count;
end

if(bottomFixed==1)
    Essential = 1:row*2;
       Essential = unique(Essential);  
end



% if(strcmp(fixedElementsCase,'bottomCorners'))
%
% elseif(strcmp(fixedElementsCase,'sidesMiddle'))%
%     Essential = []; %
%     Essential=[Essential row*2+(ceil(column/2)*row*2) 1+row*2+(ceil(column/2)*row*2)] ; % fixed at the heat source in the middle right
%     Essential=[Essential 1+(ceil(column/2)*row*2) 2+(ceil(column/2)*row*2)] ; % fixed at the heat source in the middle left
%     Essential = unique(Essential);
%     F( (floor(row/2)+1)*2) = FappliedLoad; % force down in the bottom middle
%
% elseif(strcmp(fixedElementsCase,'middleDown'))
%     Essential = []; %
%
%     Essential=[Essential row*2+(ceil(column/2)*row*2) 1+row*2+(ceil(column/2)*row*2)] ; % fixed at the heat source in the middle right
%     Essential=[Essential 1+(ceil(column/2)*row*2) 2+(ceil(column/2)*row*2)] ; % fixed at the heat source in the middle left
%     Essential = unique(Essential);
%
%     % even DOF should be y directions
%     forceNodes = [ 2+(ceil(column/2)*row*2):2:2+row*2+(ceil(column/2)*row*2)];
%     F(forceNodes)=  FappliedLoad; % middle row, downward force
%
% elseif(strcmp(fixedElementsCase,'leftClampedmiddleLoadRight'))
%

% elseif(strcmp(fixedElementsCase,'bridge'))
% %     tt=   1:2*row :2*row*(column); % ... % left side
% %     t2=tt+1;
% %     left, middle, right on bottom
% %     load down on middle row.
% %     left = [1 2];
% % %      middle = [floor(row/2)*2 floor(row/2)*2+1];
% % middle = [];
% %     right = [row*2 row*2+1];
%     nodeLeft =  floor(column/2)*row+1;
%      nodeRight =  (floor(column/2)+1)*row;
%     left =  nodeLeft*2;
%     right = nodeRight*2;
% % %     tt = leftynode:2:rightynode;
% %     right = ndof;
% %     left = ndof-row*2;
%     tt = left:2:right;
%     [~, scale] = size(tt);
%    % scale = scale(
%       Essential=[left (left -1)  right (right-1)];
%     Essential = unique(Essential);
% %     shiftupLeft = floor(column/2)*row*2+3;
% %      shiftupRight = (floor(column/2)+1)*row*2+1;
% %      tt = shiftupLeft:2:shiftupRight;
%     F(tt) = -FappliedLoad/scale;
% end



alldofs     = [1:ndof];
Free    = setdiff(alldofs,Essential);

% Add source node in the middle. Keep the temperature constant
%xMiddle  = floor(2*numNodesInRow/2);
%yMiddle = floor(2*numNodesInColumn/2);
%nodeNumber = 2*(numNodesInRow)*floor( (numNodesInColumn)/2);
% F(nodeNumber,1) = -1; % force at particular node

% F(2*(numNodesInColumn-1)*numNodesInRow+2,1) = -10; % y direction, down

%FappliedLoad = -200000;



% ke = elK_elastic(matProp);

% xLoc = 1;
% yLoc = 1;
% % loop over the elements
for e = 1:ne
    
    % loop over local node numbers to get their node global node numbers
    for j = 1:4
        % Get the node number
        coordNodeNumber = designVars.IEN(e,j);
        % get the global X,Y position of each node and put in array
        coord(j,:) = designVars.globalPosition(coordNodeNumber,:);
    end
    
    [x,y]= designVars.GivenNodeNumberGetXY(e);
    if(settings.doUseMultiElePerDV) % if elements per design var.
        [x,y] = designVars.GetDesignVarPositionGivenXYElement(settings,x,y);
    end
    
    % Get the element K matrix for this partiular element
    if(settings.macro_meso_iteration>1)
        %e = count;
        Dgiven =matProp.GetSavedDMatrix(e);
    else
        Dgiven = [];
    end
    
    [ke, KexpansionBar] = matProp.effectiveElasticKEmatrix(designVars.w(y,x), settings,Dgiven);
    % [ke] = matProp.effectiveElasticKEmatrix(  designVars.w(y,x), settings);
    
    
    
    
    
    % Insert the element stiffness matrix into the global.
    nodes1 = designVars.IEN(e,:);
    xNodes = nodes1*2-1;
    yNodes = nodes1*2;
    
    
    
    % I cannot use the union, or else the order get messed up. The order
    % is important. Same in the actual topology code when you are
    % calculating the objectiv
    NodeNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
    
    % for the x location
    % The first number is the row - "y value"
    % The second number is the column "x value"
    K(NodeNumbers,NodeNumbers) = K(NodeNumbers,NodeNumbers) + designVars.x(y,x)^settings.penal*ke;
    
    if(settings.addThermalExpansion ==1)
        alpha = matProp.effectiveThermalExpansionCoefficient(designVars.w(y,x))*designVars.x(y,x)^settings.penal;
        U_heat = designVars.U_heatColumn(nodes1,:);
        averageElementTemp = mean2(U_heat); % calculate the average temperature of the 4 nodes
        deltaTemp = averageElementTemp- settings.referenceTemperature;
        f_temperature = alpha*deltaTemp*KexpansionBar;
        F(NodeNumbers) = F(NodeNumbers) + f_temperature;
    end
    
    %        xLoc = xLoc+1;
    %        if(xLoc>settings.nelx)
    %            xLoc = 1;
    %            yLoc = yLoc+1;
    %        end
    
    
end

K = sparse(K);
F = sparse(F);
K_ff = K(Free,Free);
K_fe = K(Free,Essential);
F_f = F(Free);
% maxF = max(max(F))

% http://www.mathworks.com/help/distcomp/gpuarray.html
% http://www.mathworks.com/matlabcentral/answers/63692-matlab-cuda-slow-in-solving-matrix-vector-equation-a-x-b

if(settings.useGPU ==1)
    % GPU matrix solve.
    K_ff_gpu = gpuArray(K_ff);
    F_f_gpu = gpuArray(F_f);
    T_gpu = K_ff_gpu\F_f_gpu;
    T(Free) = gather(T_gpu);
else
    % normal matrix solve
    
    T(Free) = K_ff \ F_f;
    
end

maxF = full(max(abs(F)));
T(Essential) = u0;

maxT = full(max(T));

% plotForces =1;
% if(plotForces ==1)
%     subplot(2,2,4);
%     quiver(reshape(designVars.XLocations, ndof/2,1),reshape(designVars.YLocations, ndof/2,1),F(1:2:ndof),F(2:2:ndof));
% end

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
%  subplot(2,2,2)
%
% % % % % loop over the elements
% %  for e = 1:ne
% % %
% % %     coord = zeros(3,2);
% % %     xsum = 0;
% % %     ysum = 0;
% % %     % loop over local node numbers to get their node global node numbers
% %      for j = 1:4
% %           % Get the node number
% %           coordNodeNumber = designVars.IEN(e,j);
% %            % get the global X,Y position of each node and put in array
% %            coord(j,:) = designVars.globalPosition(coordNodeNumber,:);
% %       end
% %
% %     %  [x,y]= designVars.GivenNodeNumberGetXY(e);
% % %     elemcenterLocations(e,:) = [xsum/4 ysum/4];
% %
% %     % see version 3 of notes page 13 Also, see version 5 of the notes page
% %     % 27
% % %
% % %     eta = 0; Zeta = 0; % We are at the center, so both are zero
% % %
% % %     % B_hat (Derivative of N1 with respect to zeta and eta)
% % %      B_hat = 1/4*[-(1-eta) (1-eta) (1+eta) -(1+eta);
% % %                   -(1-Zeta) -(1+Zeta) (1+Zeta) (1-Zeta)];
% % %
% % %      % Calculate the Jacobian
% % %      J=B_hat*coord;
% % %
% % %      % Calculate the determinate
% % %      %J_det = det(J);
% % %      J_transpose = transpose(J);
% % %      J_transpose_inv = inv(J_transpose);
% % %
% % %      % Form the B matrix
% % %      B = J_transpose_inv*B_hat;
% % %
% % %     qLocal = -kmaterial*B*local_t';
% % %     qstored(e,:) = qLocal';
% % %
% % %     qMag = (qLocal(1)^2+qLocal(2)^2)^(1/2);
% % %     qMag_stored(e) =qMag;
% %
% %     % plot the element outline
% %      hold on
% %      coord(5,:) = coord(1,:);
% %      plot(coord(:,1),coord(:,2));
% %
% %  end
% % %
% %
% % % %xlabel('radial distance, meters') % y-axis label
% % % %ylabel('Height') % x-axis label
% % % tti= strcat('Flux from each element . Number of elements=', int2str(ne));
% % % title(tti);
% %   hold off
%
% q_mags = [qMag_stored,transpose(1:nn)]
%

% TcontourMatrix = zeros(numNodesInRow,numNodesInColumn);
%
% TcontourMatrixX = zeros(numNodesInRow,numNodesInColumn);
% TcontourMatrixY = zeros(numNodesInRow,numNodesInColumn);
%
%
%
% for j = 1:numNodesInColumn % y
%       rowMultiplier = j-1;
%      for i = 1:numNodesInRow % x
%          nodeNumber = i+numNodesInRow*rowMultiplier;
%          % combined the X and Y into 1 magnitude of a vector
%          node1 = nodeNumber*2-1; % minus 1, because we want to start at 1 and not 2
%          node2 = node1+1;
%          TcontourMatrixX(i,j) = T(node1);
%          TcontourMatrixY(i,j) = T(node2);
%
%          temp = sqrt(T(node1)^2+T(node2)^2);
%
%
%          TcontourMatrix(i,j) = temp;
%
%      end
% end
%
% if(1==0)
%      subplot(1,2,1)
%      % plot the coutour graph
%      contour(XLocations,YLocations,TcontourMatrix);
%      tti= strcat('Heat contours. Number of elements=', int2str(ne));
%      title(tti);
%
%     %
%     % plot the surf graph
%     subplot(1,2,2)
%     surf(XLocations,YLocations,TcontourMatrix);
% end
%  subplot(2,2,1);
%   surf(XLocations,YLocations,TcontourMatrix);
%
% if(1==1)
%      subplot(2,2,2);
%      surf(XLocations,YLocations,TcontourMatrixX);
%       subplot(2,2,3);
%      surf(XLocations,YLocations,TcontourMatrixY);
% end
% % tti= strcat('Heat surface.  Number of elements =', int2str(ne));
% % title(tti);
% %
% % diary off
T = transpose(T);
