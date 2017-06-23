function [T, maxF,maxT]=FE_elasticV2(DV, config, matProp, loadingCase)

u0 =0; % value at essentail boundaries
nn = (config.nelx+1)*(config.nely+1); % number of nodes
ne = config.nelx*config.nely; % number of elements
ndof = nn*2; % Number of degrees of freedome. 2 per node.

% Specifiy the constrained nodes where there are essential boundary
% conditions
F = zeros(ndof,1);
K = zeros(ndof,ndof);
row = config.nelx+1;
column= config.nely+1;

bottomFixed=0;

bridgeStart = (row*(column-1)+1)*2; % (row*(floor(column/2))+1)*2;
bridgeEnd=(row*(column-0))*2; % (row*(1+floor(column/2)))*2;

loadMagnitude = 10000;%500
% loading condition
if loadingCase == 111
    FappliedLoad = loadMagnitude;
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
    FappliedLoad = -loadMagnitude;
    tt=   1:2*row :2*row*(column); % ... % left side
    t2=tt+1;
    Essential=[tt t2];
    Essential = unique(Essential);    
    % down on bottom right
    F(row*2) = FappliedLoad;
    
elseif loadingCase ==113
    FappliedLoad = loadMagnitude;
    tt=   1:2*row :2*row*(column); % ... % left side
    t2=tt+1;
    Essential=[tt t2];
    Essential = unique(Essential);    
    if(mod(column,2)==0)
        % even number of column
        F( (row*2)*floor(column/2)) = FappliedLoad/2; % middle to the right  
         F( (row*2)*(floor(column/2)+1)) = FappliedLoad/2; % middle to the right  
    else
        % odd number of columns
         F( (row*2)*floor(column/2)) = FappliedLoad; % middle to the right  
    end
    

    % -------------------------------------------------
    %
    % 300s are the shoe loading case. 
    %
    % -------------------------------------------------
elseif loadingCase ==300
      FappliedLoad = -loadMagnitude;   
      bottomFixed = 1;
     
      fnodes =  (row*(column-1)+1)*2:2:(row*(column-1)+floor(row/4))*2; % top left quarter, y degrees of freedom, pointing down. 
      count = size(fnodes,2);
       F(fnodes) =FappliedLoad/count;
      
elseif loadingCase ==301
        FappliedLoad = -loadMagnitude;   
      bottomFixed = 1;
     
      fnodes =   (row*(column-1)+1)*2:2:(row*(column-1)+floor(row/2))*2; % top left half, y degrees of freedom, pointing down. 
      count = size(fnodes,2);
       F(fnodes) =FappliedLoad/count;
elseif loadingCase ==302
      FappliedLoad = -loadMagnitude;   
      bottomFixed = 1;
     
     fnodes= (row*(column-1)+1)*2:2:ndof; % whole top, y degrees of freedom, pointing down. 
      count = size(fnodes,2);
       F(fnodes) =FappliedLoad/count;
elseif loadingCase ==303
         FappliedLoad = -loadMagnitude;   
      bottomFixed = 1;
     
      fnodes =  (row*(column-1)+1+floor(row/2))*2:2:ndof; % top right half, y degrees of freedom, pointing down. 
      count = size(fnodes,2);
       F(fnodes) =FappliedLoad/count;
elseif loadingCase ==304
     FappliedLoad = -loadMagnitude;   
      bottomFixed = 1;
     
      fnodes =  (row*(column-1)+1+floor(3*row/4))*2:2:ndof; % top right quarter, x and y degrees of freedom, pointing down and left. 
      count = size(fnodes,2);
       F(fnodes) =FappliedLoad/count;
       
elseif loadingCase ==305
     FappliedLoad = -loadMagnitude;   
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
       
     FappliedLoad = -loadMagnitude;   
     bottomFixed = 1;
        
     startOffset = 0;
     endOfset = -(floor(3*row/4))*2;
        
     fnodes= bridgeStart+startOffset :2:bridgeEnd+endOfset ; % middle, 1/4 load down
     count = size(fnodes,2);
     F(fnodes) =FappliedLoad/count;

elseif loadingCase ==401
     
     FappliedLoad = -loadMagnitude;   
      bottomFixed = 1;
      
%         bridgeStart = (row*(floor(column/2))+1)*2;
%         bridgeEnd=(row*(1+floor(column/2)))*2;
        
        startOffset =(floor(1*row/4))*2;
        endOfset = -(floor(2*row/4))*2;
        
     fnodes= bridgeStart+startOffset :2:bridgeEnd+endOfset ; % middle, 1/4 load down
     count = size(fnodes,2);
     F(fnodes) =FappliedLoad/count;
elseif loadingCase ==402
      FappliedLoad = -loadMagnitude;   
      bottomFixed = 1;
      
%         bridgeStart = (row*(floor(column/2))+1)*2;
%         bridgeEnd=(row*(1+floor(column/2)))*2;
        
        startOffset =(floor(2*row/4))*2;
        endOfset = -(floor(1*row/4))*2;
        
     fnodes= bridgeStart+startOffset :2:bridgeEnd+endOfset ; % middle, 1/4 load down
     count = size(fnodes,2);
     F(fnodes) =FappliedLoad/count;
      
elseif loadingCase ==403
      FappliedLoad = -loadMagnitude;   
      bottomFixed = 1;
      
%         bridgeStart = (row*(floor(column/2))+1)*2;
%         bridgeEnd=(row*(1+floor(column/2)))*2;
        
        startOffset =(floor(3*row/4))*2;
        endOfset = 0;
        
     fnodes= bridgeStart+startOffset :2:bridgeEnd+endOfset ; % middle, 1/4 load down
     count = size(fnodes,2);
     F(fnodes) =FappliedLoad/count;
     
 elseif loadingCase ==404
      FappliedLoad = -loadMagnitude;   
      bottomFixed = 1;
      
%         bridgeStart = (row*(floor(column/2))+1)*2;
%         bridgeEnd=(row*(1+floor(column/2)))*2;
        
        startOffset =-1;
        endOfset = 0;
        
     fnodes= bridgeStart+startOffset :2:bridgeEnd+endOfset ; % middle, 1/4 load down
     count = size(fnodes,2);
     F(fnodes) =FappliedLoad/count;
  elseif loadingCase ==405
      FappliedLoad = loadMagnitude;   
      bottomFixed = 1;
      
%         bridgeStart = (row*(floor(column/2))+1)*2;
%         bridgeEnd=(row*(1+floor(column/2)))*2;
        
        startOffset =-1;
        endOfset = 0;
        
     fnodes= bridgeStart+startOffset :2:bridgeEnd+endOfset ; % middle, 1/4 load down
     count = size(fnodes,2);
     F(fnodes) =FappliedLoad/count;
end

if(bottomFixed==1)
    Essential = 1:row*2;
       Essential = unique(Essential);  
end

alldofs     = [1:ndof];
Free    = setdiff(alldofs,Essential);
E12 = 1;
E33 = 1;

% % loop over the elements
for e = 1:ne    
    % loop over local node numbers to get their node global node numbers
    for j = 1:4
        % Get the node number
        coordNodeNumber = DV.IEN(e,j);
        % get the global X,Y position of each node and put in array
        coord(j,:) = DV.globalPosition(coordNodeNumber,:);
    end
    
    [x,y]= DV.GivenNodeNumberGetXY(e);
  
    
    % Get the element K matrix for this partiular element
%     if(config.macro_meso_iteration>1)
%         %e = count;
%         Dgiven =matProp.GetSavedDMatrix(e);
%     else
%         Dgiven = [];
%     end
%     
%      if(config.useOrthAndRot==1)
%         Dgiven= matProp.calculateEffConsMatrixWithGradAndOrthDistrbution( DV.w(y,x), config, DV.d(y,x));
%     end
    
%     [ke, KexpansionBar] = matProp.effectiveElasticKEmatrix(DV.w(y,x), config,Dgiven);
    % [ke] = matProp.effectiveElasticKEmatrix(  DV.w(y,x), config);
    
    % Insert the element stiffness matrix into the global.
    nodes1 = DV.IEN(e,:);
    xNodes = nodes1*2-1;
    yNodes = nodes1*2;
    
    
    
    % I cannot use the union, or else the order get messed up. The order
    % is important. Same in the actual topology code when you are
    % calculating the objectiv
    NodeNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
    
    % for the x location
    % The first number is the row - "y value"
    % The second number is the column "x value"
%     K(NodeNumbers,NodeNumbers) = K(NodeNumbers,NodeNumbers) + DV.x(y,x)^config.penal*ke;

      % KE = matProp.getKMatrixUseTopGradOrthoDistrRotVars(config,DV.x(y,x),DV.w(y,x),DV.d(y,x),DV.t(y,x));
      if(config.anisotropicMat==1)
          E12=DV.E12(y,x);
          E33=DV.E33(y,x);
      end
   
      KE = matProp.getKMatrixTopExxYyyRotVars(config,DV.x(y,x),DV.Exx(y,x), DV.Eyy(y,x),DV.t(y,x),DV.w(y,x),E12, E33, e);
      K(NodeNumbers,NodeNumbers) = K(NodeNumbers,NodeNumbers) + KE;
    
%     if(config.addThermalExpansion ==1)
%         alpha = matProp.effectiveThermalExpansionCoefficient(DV.w(y,x))*DV.x(y,x)^config.penal;
%         U_heat = DV.U_heatColumn(nodes1,:);
%         averageElementTemp = mean2(U_heat); % calculate the average temperature of the 4 nodes
%         deltaTemp = averageElementTemp- config.referenceTemperature;
%         f_temperature = alpha*deltaTemp*KexpansionBar;
%         F(NodeNumbers) = F(NodeNumbers) + f_temperature;
%     end
    
    %        xLoc = xLoc+1;
    %        if(xLoc>config.nelx)
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

% if(config.useGPU ==1)
%     % GPU matrix solve.
%     K_ff_gpu = gpuArray(K_ff);
%     F_f_gpu = gpuArray(F_f);
%     T_gpu = K_ff_gpu\F_f_gpu;
%     T(Free) = gather(T_gpu);
% else
    % normal matrix solve
    
    T(Free) = K_ff \ F_f;
    
% end

maxF = full(max(abs(F)));
T(Essential) = u0;

maxT = full(max(T));

% plotForces =1;
% if(plotForces ==1)
%     subplot(2,2,4);
%     quiver(reshape(DV.XLocations, ndof/2,1),reshape(DV.YLocations, ndof/2,1),F(1:2:ndof),F(2:2:ndof));
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
% %           coordNodeNumber = DV.IEN(e,j);
% %            % get the global X,Y position of each node and put in array
% %            coord(j,:) = DV.globalPosition(coordNodeNumber,:);
% %       end
% %
% %     %  [x,y]= DV.GivenNodeNumberGetXY(e);
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
