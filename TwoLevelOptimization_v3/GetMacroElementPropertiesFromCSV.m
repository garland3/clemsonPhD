function  macroEleProps = GetMacroElementPropertiesFromCSV(config,e)
% for the new meso design method using the Exx and Eyy and the consistency
% constraints I need
% 1. Topology var. Tells me if I need to make a new design or not. DONE
% 2. XY position map, DONE
% 3. D_system matrix for each element in the macro design.
% 4. Displacement field. So that I can calculate the strain on each. DONE
% element.


if(config.validationModeOn==0)    
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
    
    if(macroEleProps.densitySIMP>config.noNewMesoDesignDensityCutOff)
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
    ExxValues= csvread(outname);
    
    % read the Eyy field
    outname = sprintf('./out%i/EyyValues%i.csv',folderNum,mm_iteration);
    EyyValues= csvread(outname);
    
    % read the Theta field
    outname = sprintf('./out%i/ThetaValues%i.csv',folderNum,mm_iteration);
    ThetaValues=csvread(outname);
    
    Exx=ExxValues(e);
    Eyy=EyyValues(e);
    theta=ThetaValues(e);
    
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
    
    
    
end

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