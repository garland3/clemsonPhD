function [designVars, D_h, objective,new_muMatrix]=Homgenization(designVars, config, matProp,macroElemProps,mesoLoop,old_muMatrix,penaltyValue)

% ---------------------
% Use the wrap around FEA
% -------------------
strainMultiplier = 1;

% u0 =0; % value at essentail boundaries
% nn = (config.nelx+1)*(config.nely+1); % number of nodes
nn = (config.nelx)*(config.nely); % number of nodes wrap arround.
ne = config.nelx*config.nely; % number of elements
ndof = nn*2; % Number of degrees of freedome. 2 per node.

% Specifiy the constrained nodes where there are essential boundary
% conditions
F1 = zeros(ndof,1);
F2 = zeros(ndof,1);
F3 = zeros(ndof,1);


K = zeros(ndof,ndof);
Essential = [1 2 4]; % at least 3 nodes need to be constrained

alldofs     = 1:ndof;
Free    = setdiff(alldofs,Essential);


strain1 =  [ 1 0 0]*strainMultiplier;
strain2 =  [ 0 1 0]*strainMultiplier;
strain3 =  [ 0 0 1]*strainMultiplier;

matvolFraction = 1;

% Dideal = [ 1 matProp.v 0;
%     matProp.v 1 0;
%     0 0 1/2*(1-matProp.v)]*matProp.E_material1/(1-matProp.v^2);

% Dideal=zeros(3,3);
[~, t2] = size(config.loadingCase);
% for loadcaseIndex = 1:t2
%     Dideal  = Dideal+abs(macroElemProps.strain(:,loadcaseIndex)*macroElemProps.strain(:,loadcaseIndex)'  );
% end
% Dideal=Dideal/t2;


[~, ~, B_total, ~] = matProp.effectiveElasticKEmatrix_meso(matvolFraction, config,'');


for e = 1:ne
    
    
    [x,y]= designVars.GivenNodeNumberGetXY(e);
    % The SIMP var is added later.
    [ke, ~, B_total, ~] = matProp.effectiveElasticKEmatrix_meso(matvolFraction, config,strain1);
    [~, ~, ~, F_meso1] = matProp.effectiveElasticKEmatrix_meso(matvolFraction, config,strain1);
    [~, ~, ~, F_meso2] = matProp.effectiveElasticKEmatrix_meso(matvolFraction, config,strain2);
    [~, ~, ~, F_meso3] = matProp.effectiveElasticKEmatrix_meso(matvolFraction, config,strain3);
    
    % Insert the element stiffness matrix into the global.
    nodes1 = designVars.IEN(e,:);
    xNodes = nodes1*2-1;
    yNodes = nodes1*2;
    
    % I cannot use the union, or else the order get messed up. The order
    % is important. Same in the actual topology code when you are
    % calculating the objectiv
    NodeNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
    
    
    % The constutive matrix should change based on the element's
    % topology density, so we need to apply the SIMP
    K(NodeNumbers,NodeNumbers) = K(NodeNumbers,NodeNumbers) + designVars.x(y,x)^config.penal*ke;
    F1(NodeNumbers) = F1(NodeNumbers) +F_meso1 * designVars.x(y,x)^config.penal;
    F2(NodeNumbers) = F2(NodeNumbers) +F_meso2* designVars.x(y,x)^config.penal;
    F3(NodeNumbers) = F3(NodeNumbers) +F_meso3* designVars.x(y,x)^config.penal;
    
    %
    %     if(config.addThermalExpansion ==1)
    %         alpha = matProp.effectiveThermalExpansionCoefficient(designVars.w(y,x))*designVars.x(y,x)^config.penal;
    %         U_heat = designVars.U_heatColumn(nodes1,:);
    %         averageElementTemp = mean2(U_heat); % calculate the average temperature of the 4 nodes
    %         deltaTemp = averageElementTemp- config.referenceTemperature;
    %         f_temperature = alpha*deltaTemp*KexpansionBar;
    %         F1(NodeNumbers) = F1(NodeNumbers) + f_temperature;
    %     end
    
end

K = sparse(K);
F1 = sparse(F1);
F2 = sparse(F2);
F3 = sparse(F3);
F_f1 = F1(Free);
F_f2 = F2(Free);
F_f3 = F3(Free);
K_ff = K(Free,Free);
% K_fe = K(Free,Essential);
% http://www.mathworks.com/help/distcomp/gpuarray.html
% http://www.mathworks.com/matlabcentral/answers/63692-matlab-cuda-slow-in-solving-matrix-vector-equation-a-x-b

if(config.useGPU ==1)
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

u0=0;
T1(Essential) = u0;
T2(Essential) = u0;
T3(Essential) = u0;

% D constriutive matrix, homoegenized, but the sum, not averaged yet.
D_h = zeros(3,3);

designVars.dc = zeros(config.nely,config.nelx);


objective = 0;

muVersion2 = ones(1,3);
%muTarget =config. muSequence(  config.macro_meso_iteration);
term1Max = 0;
term2Max = 0;

%    if(config.addConsistencyConstraints==1 &&mesoLoop>2 )
Diff_Sys_Sub =  (macroElemProps.D_subSys- macroElemProps.D_sys);

% Diff_Sys_Sub = Diff_Sys_Sub/max(max(Diff_Sys_Sub)); % scale down.
%    end
% dx3 = zeros(config.nely,config.nelx);

% targetStrain = inv( macroElemProps.D_sys)*ones(3,1)*10;

for e = 1:ne
    [x,y]= designVars.GivenNodeNumberGetXY(e);
    nodes1=  designVars.IEN(e,:);
    xNodes = nodes1*2-1;
    yNodes = nodes1*2;
    dofNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
    Ulocal1 = T1(dofNumbers);
    Ulocal2 = T2(dofNumbers);
    Ulocal3 = T3(dofNumbers);
    %         material1Fraction = designVars.w(y,x); % 100% of material 1 right now.
    material1Fraction=1;
    E_base =    matProp.effectiveElasticProperties( material1Fraction, config);
    %     E = E_base;
    
    v = 0.3; % Piossons ratio
    % D is called C* in some journal papers.
    D_base = [ 1 v 0;
        v 1 0;
        0 0 1/2*(1-v)]*E_base/(1-v^2);
    
    Ulocal1 = full(Ulocal1);Ulocal2 = full(Ulocal2);Ulocal3 = full(Ulocal3);
    temp1_X = [transpose(Ulocal1) transpose(Ulocal2) transpose(Ulocal3)];
    temp2_BX = B_total*temp1_X;
    temp3 = (eye(3)*strainMultiplier-temp2_BX);
    
    D = D_base*designVars.x(y,x)^config.penal; % add the simp multipler
    D_h_element = transpose(temp3)*D*temp3;
    %       D_h_element = D*temp3;
    D_h = D_h_element+D_h;
    
    A = 3/ne*transpose(temp3)*D_base*temp3;
    
    
    % Calculate the sensitivity
    D_d = D_base*config.penal*  designVars.x(y,x)^(config.penal-1);
    dH =  transpose(temp3)*D_d*temp3;
    %     dH =  D_d*temp3;
    
    %      D_d3 = D_base*config.penal*(config.penal-1)*  designVars.x(y,x)^(config.penal-2);
    %      dH3= transpose(temp3)*D_d3*temp3;
    
    
    % inverse homogenization maximize the shear stiffness
    % designVars.dc(y,x) =-dH(3,3);
    
    % maxmize the bulk modulus
    % designVars.dc(y,x) =-(dH(1,1)+dH(2,2)+ dH(1,2)+dH(2,1)); %
    
    % maximize the stiffness in the y direciton
    % designVars.dc(y,x) = -dH(2,2);
    mode =5;
    
    % ----------------------------------
    % Update the sensitiviy
    % mode 1 = playing around
    % mode 2 = mu term to inforce consistency (not working), the second
    % term can only contribute a fraction amount
    % mode 3 = mu and lagrangian term using augmented lagragian
    % mode 4 = mu and lagrangian terms where using a time step
    % mode 5 = find equavilent epsilon (strain), then maximize strain
    % energy
    % ----------------------------------
    if(mode ==1)
        designVars.dc(y,x) =-dH(1,1);
        
    elseif(mode==2)
        % two scale optimization
        for loadcaseIndex = 1:t2
            designVars.dc(y,x)  = -macroElemProps.strain(:,loadcaseIndex)'*dH*macroElemProps.strain(:,loadcaseIndex)+ designVars.dc(y,x)  ;
        end
        term1Max = max(term1Max  ,abs(designVars.dc(y,x)));
        % This might should be added in the optimal criteria section, but
        % I'm trying it here first.
        if(config.addConsistencyConstraints==1 &&mesoLoop>2 )
            term2 = 2*muVersion2*(Diff_Sys_Sub.*dH.*sign(Diff_Sys_Sub))*transpose(muVersion2);
            term2Max=max(term2Max,abs(term2));
            term2=muTarget*term2*oldTerm1Max/oldTerm2Max;
            designVars.dc(y,x)=  designVars.dc(y,x)+term2;
        else
            term2Max=1000;
        end
    elseif(mode ==3)
        % Auguemented Lagrangian Multiplier optimization
        % Solve Equation 19
        rowIndex = [1,1,2,3];
        columnIndex = [1,2,2,3];
        
        P  =0;
        %       Ptemp=0;
        constraintCount = 0;
        for k = 1:4
            i = rowIndex(k);
            j = columnIndex(k);
            % temp = old_muMatrix(i,j)*-dH(i,j)+2*penaltyValue*Diff_Sys_Sub(i,j)*-dH(i,j);
            Ptemp = A(i,j)*(-old_muMatrix(i,j)-penaltyValue*Diff_Sys_Sub(i,j));
            P =P +Ptemp;
            constraintCount=constraintCount+1;
        end
        move = 0.1;
        optimalEta = -2/(constraintCount*P);
        xx=designVars.x(y,x); % =min(optimalEta, designVars.x+move)
        
        designVars.x(y,x) = max(0.01,max(xx-move,min(1.,min(xx+move,optimalEta))));
        
        
        
        
        designVars.dc(y,x) =dH(1,1);
        
    elseif(mode ==4)
       
           designVars.d11(y,x) =dH(1,1);
           designVars.d12(y,x) =dH(1,2);
           designVars.d22(y,x) =dH(2,2);
           designVars.d33(y,x) =dH(3,3);
           
           designVars.De11(y,x) =D_d(1,1);
           designVars.De12(y,x) =D_d(1,2);
           designVars.De22(y,x) =D_d(2,2);
           designVars.De33(y,x) =D_d(1,3);
           
            designVars.dc(y,x) =dH(1,1)*old_muMatrix(1,1) +  dH(1,2)*old_muMatrix(1,2) ;% +  dH(2,2)*old_muMatrix(2,2)  +dH(3,3)*old_muMatrix(3,3);
    elseif(mode ==5)
%         loadcaseIndex=1;
          designVars.dc(y,x)  = -macroElemProps.psuedoStrain'*dH*macroElemProps.psuedoStrain;
%           tempStrain = macroElemProps.psuedoStrain;
%           tempStrain(3)=tempStrain(3)*-1;
%            designVars.dc(y,x)  = -tempStrain'*dH*tempStrain+ designVars.dc(y,x) ;
        
    end
    
end


%lambdaMatrix = lambdaMatrix-muTarget*(Dideal- macroElemProps.D_sys);

% lambdaMatrix = lambdaMatrix-muTarget*(D_h- macroElemProps.D_sys);
% sys =  macroElemProps.D_sys;
D_h = D_h/ne;
% tempDiff =sys- D_h;
% damping = 0.05;
% new_muMatrix=old_muMatrix+damping*penaltyValue*(tempDiff);
new_muMatrix=old_muMatrix;

% now that we have D_h, calculate the energy objective.
for loadcaseIndex = 1:t2
    loadcaseIndex=1;
    objective =objective+ macroElemProps.strain(:,loadcaseIndex)'*D_h*macroElemProps.strain(:,loadcaseIndex);
end

% force to be negative.
% maxdc = max(max( designVars.dc));
% if(maxdc>-0.001)
%     designVars.dc =  designVars.dc-maxdc-1;
% end
%
% % auto scale to be  between -1 and -100
%
%  largest = max(max(designVars.dc));
%  smallest = min(min(   designVars.dc));
%     designVars.dc=    designVars.dc-largest;
%      designVars.dc= designVars.dc*100/(abs(smallest-largest));
%designVars.dc =  designVars.dc-50000;
