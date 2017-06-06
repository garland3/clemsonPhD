function  macroEleProps = GetMacroElementPropertiesFromCSV(config,e)
% for the new meso design method using the Exx and Eyy and the consistency
% constraints I need
% 1. Topology var. Tells me if I need to make a new design or not. DONE
% 2. XY position map, DONE
% 3. D_system matrix for each element in the macro design.
% 4. Displacement field. So that I can calculate the strain on each. DONE
% element.


if(config.validationModeOn==0 && config.strainAndTargetTest==0)
    % ----------------------------------------------
    %
    %      Normal Multiscale optimization
    %
    % ----------------------------------------------
    
    mm_iteration = config.macro_meso_iteration;
    macroEleProps = macroElementProp;
    macroEleProps.elementNumber = e;
    
    folderNum = config.iterationNum;
    
    % Get element->node mapping
    outname = sprintf('./out%i/elementNodeMap%i.csv',folderNum,mm_iteration);
    IEN = csvread(outname);
    
    % % Get displacement fields
    outname = sprintf('./out%i/displacement%i.csv',folderNum,mm_iteration);
    U =  csvread(outname);
    
    % -----------------------------------
    %
    % 1 element per design var case on macro level
    %
    % -----------------------------------
    % GET the saved element to XY position map (needed for x and w vars retrival)
    outname = sprintf('./out%i/elementXYposition%i.csv',folderNum,mm_iteration);
    elementXYposition=csvread(outname);
    results = elementXYposition(macroEleProps.elementNumber,:);
    macroEleProps.yPos = results(1);
    macroEleProps.xPos = results(2);
    
    % Get the density field
    outname = sprintf('./out%i/SIMPdensityfield%i.csv',folderNum,mm_iteration);
    x = csvread(outname);
    macroEleProps.densitySIMP = x(macroEleProps.yPos,macroEleProps.xPos );
    
    % Get the Exx field
    outname = sprintf('./out%i/ExxValues%i.csv',folderNum,mm_iteration);
    ExxMacro = csvread(outname);
    macroEleProps.Exx = ExxMacro(macroEleProps.yPos,macroEleProps.xPos );
    
    % Get the Eyy field
    outname = sprintf('./out%i/EyyValues%i.csv',folderNum,mm_iteration);
    EyyMacro =csvread(outname);
    macroEleProps.Eyy = EyyMacro(macroEleProps.yPos,macroEleProps.xPos );
    
    % Get the Theta field
    outname = sprintf('./out%i/ThetaValues%i.csv',folderNum,mm_iteration);
    ThetaMacro = csvread(outname);
    macroEleProps.theta = ThetaMacro(macroEleProps.yPos,macroEleProps.xPos );
    
    w=1;
    matProp=MaterialProperties;
    D= matProp.getDmatMatrixTopExxYyyRotVars(config,macroEleProps.densitySIMP ,macroEleProps.Exx, macroEleProps.Eyy,macroEleProps.theta ,w);
    
    % outname = sprintf('./out%i/DsystemIter%i_Element_%i.csv',folderNum,mm_iteration,e);
    % D = csvread(outname);
    macroEleProps.D_sys =D;
    
    if(macroEleProps.densitySIMP>config.noNewMesoDesignDensityCutOff ||config.multiscaleMethodCompare==1)
        nodes1=  IEN(e,:);
        macroEleProps.elementNodes=nodes1;
        xNodes = nodes1*2-1;
        yNodes = nodes1*2;
        dofNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
        
        % plan for multi-loading cases.
        [~, t2] = size(config.loadingCase);
        for loadcaseIndex = 1:t2
            utemp = U(loadcaseIndex,:);
            u_local =   utemp(dofNumbers);
            
            %         offsetX = mean(u_local([1 3 5 7]));
            %         offsetY = mean(u_local([2 4 6 8]));
            offsetX = u_local(7);
            offsetY = u_local(8);
            u_local([1 3 5 7]) = u_local([1 3 5 7])-offsetX;
            u_local([2 4 6 8]) = u_local([2 4 6 8])-offsetY;
            macroEleProps.disp(loadcaseIndex,:)  = u_local;
        end
        
        % % Get the volume fraction field (but only if using the old method)
        if(config.useExxEyy~=1)
            outname = sprintf('./out%i/volfractionfield%i.csv',folderNum,mm_iteration);
            w = csvread(outname);
            macroEleProps.material1Fraction = w(macroEleProps.yPos,macroEleProps.xPos );
        else
            macroEleProps.material1Fraction = 1; % set to 1 as a place holder
        end
    else
        fprintf('SIMP density is below threshold\n');
    end
elseif(config.validationModeOn==1)
    % -----------------------------------
    %
    %       Validation problem for meso designs.
    %
    % -----------------------------------
    folderNum=0;
    mm_iteration=1;
    % read the Exx field
    outname = sprintf('./out%i/ExxValues%i.csv',folderNum,mm_iteration);
    strain1= csvread(outname);
    
    % read the Eyy field
    outname = sprintf('./out%i/EyyValues%i.csv',folderNum,mm_iteration);
    strain2= csvread(outname);
    
    % read the Theta field
    outname = sprintf('./out%i/ThetaValues%i.csv',folderNum,mm_iteration);
    ThetaValues=csvread(outname);
    
    Exx=strain1(e);
    Eyy=strain2(e);
    theta=ThetaValues(e);
    fprintf('Meso design Validation mode. Targets Exx %f Eyy %f Theta %f\n',Exx,Eyy,theta);
    
    macroEleProps = macroElementProp;
    macroEleProps.elementNumber = e;
    macroEleProps.yPos =1;
    macroEleProps.xPos = e;
    macroEleProps.material1Fraction=1;
    macroEleProps.disp=ones(1,8);
    
    
    macroEleProps.densitySIMP = 1;
    macroEleProps.Exx =Exx;
    macroEleProps.Eyy =Eyy;
    macroEleProps.theta=theta;
    
    
    w=1;
    matProp=MaterialProperties;
    D= matProp.getDmatMatrixTopExxYyyRotVars(config,macroEleProps.densitySIMP ,macroEleProps.Exx, macroEleProps.Eyy,macroEleProps.theta ,w);
    
    % outname = sprintf('./out%i/DsystemIter%i_Element_%i.csv',folderNum,mm_iteration,e);
    % D = csvread(outname);
    macroEleProps.D_sys =D;
    
    
    
elseif(config.strainAndTargetTest==1)
    % -----------------------------------
    %
    %      Psuedo Strain and target density test for ANN training data
    %
    % -----------------------------------
    folderNum=0;
    mm_iteration=1;
    % pseudoStrain
    outname = sprintf('./out%i/strain1%i.csv',folderNum,mm_iteration);
    strain1= csvread(outname);
    
    % pseudoStrain
    outname = sprintf('./out%i/strain2%i.csv',folderNum,mm_iteration);
    strain2= csvread(outname);
    
    % pseudoStrain
    outname = sprintf('./out%i/strain3%i.csv',folderNum,mm_iteration);
    strain3=csvread(outname);
    
    % targetDensity
    outname = sprintf('./out%i/densityTargets%i.csv',folderNum,mm_iteration);
    densityTarget=csvread(outname);
    
    %     Exx=strain1(e);
    %     Eyy=strain2(e);
    %     theta=ThetaValues(e);
    
    macroEleProps = macroElementProp;
    macroEleProps.elementNumber = e;
    macroEleProps.yPos =1;
    macroEleProps.xPos = e;
    macroEleProps.material1Fraction=1;
    macroEleProps.disp=ones(1,8);
    
    macroEleProps.psuedoStrain=[1; 1; 1];
    macroEleProps.psuedoStrain(1) = strain1(e);
    macroEleProps.psuedoStrain(2) = strain2(e);
    macroEleProps.psuedoStrain(3) = strain3(e);
    temp1 =densityTarget(e);
    macroEleProps.targetDensity=temp1;
    
    
    macroEleProps.densitySIMP = 1;
    macroEleProps.Exx =1;
    macroEleProps.Eyy =1;
    macroEleProps.theta=1;
    
    
    w=1;
    matProp=MaterialProperties;
    D= matProp.getDmatMatrixTopExxYyyRotVars(config,macroEleProps.densitySIMP ,macroEleProps.Exx, macroEleProps.Eyy,macroEleProps.theta ,w);
    
    % outname = sprintf('./out%i/DsystemIter%i_Element_%i.csv',folderNum,mm_iteration,e);
    % D = csvread(outname);
    macroEleProps.D_sys =D;
    
    %printf('Finished getting macro data from csv files');
end

if(config.UseLookUpTableForPsuedoStrain==1 && config.strainAndTargetTest~=1)
    
    % 1. REad in the data look up table.
    % 2. (Done) find the D_h for the coorospnding Exx, Eyy, Theta values
    % 3. Search the lookup table for the value that minimizes
    % Save this min value as the psuedo strain to use.
    macro_meso_iteration=1;
    outname = sprintf('./data/D11%i_datat.data',macro_meso_iteration);
    D11 =csvread(outname);
    
    %     outname = sprintf('./data/D12%i_datat.data',macro_meso_iteration);
    %      D12 =csvread(outname);
    
    outname = sprintf('./data/D22%i_datat.data',macro_meso_iteration);
    D22 =csvread(outname);
    
    outname = sprintf('./data/D33%i_datat.data',macro_meso_iteration);
    D33 =csvread(outname);
    
    % write the targets
    outname = sprintf('./data/pstrain1%i_datat.data',macro_meso_iteration);
    pstrain1 =csvread(outname);
    
    outname = sprintf('./data/pstrain2%i_datat.data',macro_meso_iteration);
    pstrain2 =csvread(outname);
    
    outname = sprintf('./data/pstrain3%i_datat.data',macro_meso_iteration);
    pstrain3 =csvread(outname);
    
    outname = sprintf('./data/etaTarget%i_datat.data',macro_meso_iteration);
    etaTarget =csvread(outname);
    
    
    [t1 t2]=size(etaTarget);
    
    minValue=1e11;
    indexOfMinValue = 1;
    
    D11sys = macroEleProps.D_sys(1,1) ;
    D22sys = macroEleProps.D_sys(2,2) ;
    D33sys = macroEleProps.D_sys(3,3) ;
    for i = 1:t2
        D11_table = D11(i);
        D22_table = D22(i);
        D33_table = D33(i);
        
          diffValue = (D11_table-D11sys)^2+(D22_table-D22sys)^2+(D33_table-D33sys)^2;
%           diffValue = (D11_table-D11sys)^2+(D22_table-D22sys)^2+(D33_table-D33sys)^2;
%               diffValue = (D11_table-D11sys)^2+(D22_table-D22sys)^2;
        
        if(diffValue<minValue)
            minValue=diffValue;
            indexOfMinValue=i;
        end
    end
    
    D11_table = D11(indexOfMinValue);
    D22_table = D22(indexOfMinValue);
    D33_table = D33(indexOfMinValue);
    
    diffValue = (D11_table-D11sys)^2+(D22_table-D22sys)^2;%+0.5*(D33_table-D33sys)^2;
    
    ps(1) = pstrain1(indexOfMinValue);
   ps(2)= pstrain2(indexOfMinValue);
   ps(3) = pstrain3(indexOfMinValue);
   etaTargetLocal = etaTarget(indexOfMinValue);
    
    
     macroEleProps.psuedoStrain(1) = ps(1) ;
    macroEleProps.psuedoStrain(2) = ps(2);
    macroEleProps.psuedoStrain(3) =  ps(3);
    macroEleProps.targetDensity=etaTargetLocal;
    
    
    
    
end
% else
%
%     % -----------------------------------
%     %
%     % More than 1 element per design var case on macro level
%     % e = design var and not element for this case.
%     %
%     % -----------------------------------
%
%
%     %
%     %1. make an array of meso local X positions , and make another one for Y positions for each node that is a part of the elements
%     %2. make an zeros array that will hold the X displacements at each
%     %3. make an zeros array that will hold the Y displacements at each
%     %4. make a list of macro element node numbers
%     %5. Loop over the macro elements
%     %	1. check if this node has already been added to the lsit of node numbers
%     %	2. if not then add the X and Y displacements to the local displacement array
%
%     %     num Xnodes = config.numXElmPerDV
%     %      [X,Y] = ndgrid(0:config.numXElmPerDV,0:config.numYElmPerDV);
%     [Y,X] = ndgrid(0:config.numYElmPerDV,0:config.numXElmPerDV);
%     %     [t1, t2] = size(X);
%     %     displacementsX = zeros(t1,t2);
%     %     displacementsY = displacementsX;
%
%     macroEleProps.mesoXnodelocations=X;
%     macroEleProps.mesoYnodelocations=Y;
%
%     numVarsinRow =config.numVarsX;
%     %     numVarsinColumn =config.numVarsY;
%     ydesignVar =   floor( (e-1)/numVarsinRow)+1;
%     xdesignVar = mod(e-1,numVarsinRow)+1;
%
%     % trying to find the nodes and then dofs of the macro problem
%     shiftElementsX = (xdesignVar-1)*config.numXElmPerDV;
%     elementsInRow = config.nelx;
%     shiftNodesY = (ydesignVar-1)*config.numYElmPerDV*elementsInRow;
%
%     %     count1 = 1;
%     %     count2 = 1;
%     %     isfirstrow = 1;
%
%     nodelist=[];
%     for i  =1: config.numYElmPerDV
%         startElement = shiftElementsX+shiftNodesY+(i-1)*elementsInRow;
%         for j =1: config.numXElmPerDV
%             ee = startElement+j;
%
%             nodes1=  IEN(ee,:);
%             nodelist = [nodelist nodes1];
%         end
%     end
%
%     nodelist = unique(sort(nodelist));
%
%     % nodes1=  IEN(e,:);
%     macroEleProps.elementNodes=nodelist;
%
%     % multiloading cases.
%     [~, t2] = size(config.loadingCase);
%     for loadcaseIndex = 1:t2
%         utemp = U(loadcaseIndex,:);
%
%         displacementsX = utemp(nodelist*2-1);
%         displacementsY = utemp(nodelist*2);
%
%         macroEleProps.xDisplacements(loadcaseIndex,:)=displacementsX;
%         macroEleProps.yDisplacements(loadcaseIndex,:)=displacementsY;
%     end
%
%     macroEleProps.elementNumber = e;
%     % Save element to XY position map (needed for x and w vars retrival)
%     %         outname = sprintf('./out%i/elementXYposition%i.csv',folderNum,macro_meso_iteration);
%     %         elementXYposition=csvread(outname);
%     %         results = elementXYposition(macroEleProps.elementNumber,:);
%     macroEleProps.yPos = ydesignVar;
%     macroEleProps.xPos = xdesignVar;
%
%     % Get the density field
%     outname = sprintf('./out%i/densityfield%i.csv',folderNum,macro_meso_iteration);
%     x = csvread(outname);
%     macroEleProps.density = x(macroEleProps.yPos,macroEleProps.xPos );
%
%     % % Get the volume fraction field
%     outname = sprintf('./out%i/volfractionfield%i.csv',folderNum,macro_meso_iteration);
%     w = csvread(outname);
%     macroEleProps.material1Fraction = w(macroEleProps.yPos,macroEleProps.xPos );
%
%     % if(macro_meso_iteration>1)
%     outname = sprintf('./out%i/Dgiven_%i_forElement_%i.csv',folderNum,macro_meso_iteration,e);
%     D = csvread(outname);
%     macroEleProps.D_given =D;
%     % end
%
%
%
% end