function [D_homog]=FE_elasticV2_homgonization(designVars, settings, matProp)

    


% E = matProp.E_material1; % Young's mod
% v = matProp.v; % Piossons ratio
% G = matProp.G;

u0 =0; % value at essentail boundaries
nn = (settings.nelx)*(settings.nely); % number of nodes
ne = settings.nelx*settings.nely; % number of elements
ndof = nn*2; % Number of degrees of freedome. 2 per node. 

% ------------------------------------
% -- Make a mapping between elements and global matrixes
% ------------------------------------

% IEN holds the node numbers for each element. 
% Each row is a new element
% The first column is element 1's global node number
% Second column is elemnt 2's global node number
%  and ....

%  count = 1;
%  elementsInRow = settings.nelx+1;
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
%  count = 1;
%  for i = 1:numNodesInColumn  % y
%      for j= 1:numNodesInRow % x
%          globalPosition(count,:) = [j-1 i-1];
%          count = count+1;        
%          XLocations(j,i) = j-1;
%          YLocations(j,i) = i-1;
%      end
%  end 

% Specifiy the constrained nodes where there are essential boundary
% conditions
F1 = zeros(ndof,1);
F2 = zeros(ndof,1);
F3 = zeros(ndof,1);
K = zeros(ndof,ndof);
row = settings.nelx;
column= settings.nely;
Essential = [];
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
         
%  Essential   =union( 1:2*numNodesInRow :2*(numNodesInColumn-1)*numNodesInRow+1,numNodesInRow*2) ;% X direction on left side is held fixed
               % Y direction on left side is held fixed
% Essential = unique(Essential);

% 
% 
% if(strcmp(fixedElementsCase,'bottomCorners'))
%     Essential = [1 2]; % bottom left corner is fixed
%     Essential = [Essential row*2 row*2-1]; % bottom right corner is fixed
%     Essential = unique(Essential);
%     F( (floor(row/2)+1)*2) = FappliedLoad; % force down in the bottom middle
%     
% elseif(strcmp(fixedElementsCase,'sidesMiddle'))
%     
%     Essential = []; % 
%     Essential=[Essential row*2+(ceil(column/2)*row*2) 1+row*2+(ceil(column/2)*row*2)] ; % fixed at the heat source in the middle right
%     Essential=[Essential 1+(ceil(column/2)*row*2) 2+(ceil(column/2)*row*2)] ; % fixed at the heat source in the middle left
%     Essential = unique(Essential);
%     F( (floor(row/2)+1)*2) = FappliedLoad; % force down in the bottom middle
% 
% elseif(strcmp(fixedElementsCase,'middleDown'))    
%      Essential = []; % 
%      
%      Essential=[Essential row*2+(ceil(column/2)*row*2) 1+row*2+(ceil(column/2)*row*2)] ; % fixed at the heat source in the middle right
%     Essential=[Essential 1+(ceil(column/2)*row*2) 2+(ceil(column/2)*row*2)] ; % fixed at the heat source in the middle left
%     Essential = unique(Essential);
%      
%     % even DOF should be y directions
%     forceNodes = [ 2+(ceil(column/2)*row*2):2:2+row*2+(ceil(column/2)*row*2)];
%      F(forceNodes)=  FappliedLoad; % middle row, downward force
%    
% elseif(strcmp(fixedElementsCase,'leftClamped'))   
% %      Essential = [ 1 2 ]; % 
% %      for i = 1:column-1
% %          Essential = [Essential (i*row*2)+3 (i*row*2)+4 ]; % 
% %      end
%     tt=   1:2*row :2*row*(column-1); % ... % left side
%     t2=tt+1;
%      Essential=[tt t2];
%        Essential = unique(Essential);
%        if(isinteger(column/2))
%             F((column/2)*row*2) = FappliedLoad;
%        else
%              F(ceil(column/2)*row*2) = FappliedLoad/2.0;
%              F(floor(column/2)*row*2) = FappliedLoad/2.0;
%        end
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
%strain =
strain1 =  [ 1 0 0];
strain2 =  [ 0 1 0];
strain3 =  [ 0 0 1];
  
% xLoc = 1;
% yLoc = 1;
B_total = [];
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
    
       [ke, KexpansionBar, B_total, F_meso1] = matProp.effectiveElasticKEmatrix_meso(designVars.w(y,x), settings,strain1);
     
       [~, ~, ~, F_meso2] = matProp.effectiveElasticKEmatrix_meso(designVars.w(y,x), settings,strain2);
       [~, ~, ~, F_meso3] = matProp.effectiveElasticKEmatrix_meso(designVars.w(y,x), settings,strain3);
       
      
 
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
      
       % The constutive matrix should change based on the element's
       % topology density, so we need to apply the SIMP 
       K(NodeNumbers,NodeNumbers) = K(NodeNumbers,NodeNumbers) + designVars.x(y,x)^settings.penal*ke;
       F1(NodeNumbers) = F1(NodeNumbers) +F_meso1* designVars.x(y,x)^settings.penal;
       F2(NodeNumbers) = F2(NodeNumbers) +F_meso2* designVars.x(y,x)^settings.penal;
       F3(NodeNumbers) = F3(NodeNumbers) +F_meso3* designVars.x(y,x)^settings.penal;
       
       
       if(settings.addThermalExpansion ==1)
            alpha = matProp.effectiveThermalExpansionCoefficient(designVars.w(y,x))*designVars.x(y,x)^settings.penal;
            U_heat = designVars.U_heatColumn(nodes1,:);
            averageElementTemp = mean2(U_heat); % calculate the average temperature of the 4 nodes
            deltaTemp = averageElementTemp- settings.referenceTemperature;
            f_temperature = alpha*deltaTemp*KexpansionBar;
            F1(NodeNumbers) = F1(NodeNumbers) + f_temperature;
       end
       
%        xLoc = xLoc+1;
%        if(xLoc>settings.nelx)
%            xLoc = 1;
%            yLoc = yLoc+1;
%        end
       
       
end
  
K = sparse(K);
F1 = sparse(F1);
F2 = sparse(F2);
F3 = sparse(F3);


F_f1 = F1(Free);
F_f2 = F2(Free);
F_f3 = F3(Free);

K_ff = K(Free,Free);
K_fe = K(Free,Essential);
% maxF = max(max(F))
  
% http://www.mathworks.com/help/distcomp/gpuarray.html
% http://www.mathworks.com/matlabcentral/answers/63692-matlab-cuda-slow-in-solving-matrix-vector-equation-a-x-b

if(settings.useGPU ==1)
    % GPU matrix solve. 
    K_ff_gpu = gpuArray(K_ff);
    F_f_gpu = gpuArray(F_f1);
    T_gpu = K_ff_gpu\F_f_gpu;
    T1(Free) = gather(T_gpu);
else
    % normal matrix solve
     
     T1(Free) = K_ff \ F_f1;
     T2(Free) = K_ff \ F_f2;
     T3(Free) = K_ff \ F_f3;

end

maxF = 0;
T1(Essential) = u0;
T2(Essential) = u0;
T3(Essential) = u0;

maxT = 0;

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
% % % loop over the elements

% D constriutive matrix, homoegenized, but the sum, not averaged yet. 
D_sum_h = zeros(3,3);

for e = 1:ne

          [x,y]= designVars.GivenNodeNumberGetXY(e);
          [~, ~, B_total, ~] = matProp.effectiveElasticKEmatrix_meso(designVars.w(y,x), settings,'');
     
        nodes1=  designVars.IEN(e,:);
                   % nodes1=[rowMultiplier*elementsInRow+elx;
                    %    rowMultiplier*elementsInRow+elx+1;
                     %   (rowMultiplier +1)*elementsInRow+elx+1;
                      %  (rowMultiplier +1)*elementsInRow+elx];
                    
    xNodes = nodes1*2-1;
    yNodes = nodes1*2;
    dofNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];

    Ulocal1 = T1(dofNumbers);
     Ulocal2 = T2(dofNumbers);
      Ulocal3 = T3(dofNumbers);
      
      material1Fraction = designVars.w(y,x); % 100% of material 1 right now. 
      
      E_base =    matProp.effectiveElasticProperties( material1Fraction, settings);
      
      
      E = E_base*designVars.x(y,x)^settings.penal;
          v = 0.3; % Piossons ratio
      
      % D is called C* in some journal papers. 
        D = [ 1 v 0;
        v 1 0;
        0 0 1/2*(1-v)]*E/(1-v^2);
    
    %temp1 = [Ulocal1; Ulocal2; Ulocal3]
    Ulocal1 = full(Ulocal1);Ulocal2 = full(Ulocal2);Ulocal3 = full(Ulocal3);
    temp1 = [transpose(Ulocal1) transpose(Ulocal2) transpose(Ulocal3)];
    %temp1 = full(temp1);
    temp2 = B_total*temp1;
    temp3 = (eye(3)-temp2);
    temp4 = transpose(temp3)*D*temp3;
    D_sum_h = D_sum_h+temp4;
    
 
    
end

D_h = D_sum_h/ne;
D_homog = D_h;
% T1 = transpose(T1);
