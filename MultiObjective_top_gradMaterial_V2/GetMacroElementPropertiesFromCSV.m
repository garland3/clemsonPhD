function  macroElementProps = GetMacroElementPropertiesFromCSV(settings,e)

macro_meso_iteration = settings.macro_meso_iteration;
macroElementProps = macroElementProp;
macroElementProps.elementNumber = e;

folderNum = settings.iterationNum;

% GEt element->node mapping
outname = sprintf('./out%i/elementNodeMap%i.csv',folderNum,macro_meso_iteration);IEN = csvread(outname);

% % Get displacement field
outname = sprintf('./out%i/displacement%i.csv',folderNum,macro_meso_iteration);

U =  csvread(outname);




% -----------------------------------
%
% 1 element per design var case on macro level
%
% -----------------------------------
if(settings.doUseMultiElePerDV~=1) % if elements per design var.
    
    nodes1=  IEN(e,:);
    macroElementProps.elementNodes=nodes1;
    xNodes = nodes1*2-1;
    yNodes = nodes1*2;
    dofNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
    
    % plan for multi-loading cases. 
    [~, t2] = size(settings.loadingCase);        
    for loadcaseIndex = 1:t2
        utemp = U(loadcaseIndex,:);    
        u_local =   utemp(dofNumbers);
   
        %         offsetX = mean(u_local([1 3 5 7]));
        %         offsetY = mean(u_local([2 4 6 8]));
        offsetX = u_local(7);
        offsetY = u_local(8);
        u_local([1 3 5 7]) = u_local([1 3 5 7])-offsetX;
        u_local([2 4 6 8]) = u_local([2 4 6 8])-offsetY;
        macroElementProps.disp(loadcaseIndex,:)  = u_local;
    end
    
%     macroElementProps.elementNumber
    % Save element to XY position map (needed for x and w vars retrival)
    outname = sprintf('./out%i/elementXYposition%i.csv',folderNum,macro_meso_iteration);
    elementXYposition=csvread(outname);
    results = elementXYposition(macroElementProps.elementNumber,:);
    macroElementProps.yPosition = results(1);
    macroElementProps.xPosition = results(2);
    
    % Get the density field
    outname = sprintf('./out%i/densityfield%i.csv',folderNum,macro_meso_iteration);
    x = csvread(outname);
    macroElementProps.density = x(macroElementProps.yPosition,macroElementProps.xPosition );
    
    % % Get the volume fraction field
    outname = sprintf('./out%i/volfractionfield%i.csv',folderNum,macro_meso_iteration);
    w = csvread(outname);
    macroElementProps.material1Fraction = w(macroElementProps.yPosition,macroElementProps.xPosition );
    
    % if(macro_meso_iteration>1)
    outname = sprintf('./out%i/Dgiven_%i_forElement_%i.csv',folderNum,macro_meso_iteration,e);
    D = csvread(outname);
    macroElementProps.D_given =D;
    % end
    
else
    
    % -----------------------------------
    %
    % More than 1 element per design var case on macro level
    % e = design var and not element for this case.
    %
    % -----------------------------------
    
    
    %
    %1. make an array of meso local X positions , and make another one for Y positions for each node that is a part of the elements
    %2. make an zeros array that will hold the X displacements at each
    %3. make an zeros array that will hold the Y displacements at each
    %4. make a list of macro element node numbers
    %5. Loop over the macro elements
    %	1. check if this node has already been added to the lsit of node numbers
    %	2. if not then add the X and Y displacements to the local displacement array
    
    %     num Xnodes = settings.numXElmPerDV
    %      [X,Y] = ndgrid(0:settings.numXElmPerDV,0:settings.numYElmPerDV);
    [Y,X] = ndgrid(0:settings.numYElmPerDV,0:settings.numXElmPerDV);
%     [t1, t2] = size(X);
%     displacementsX = zeros(t1,t2);
%     displacementsY = displacementsX;
    
    macroElementProps.mesoXnodelocations=X;
    macroElementProps.mesoYnodelocations=Y;
    
   
    numVarsinRow =settings.numVarsX;
%     numVarsinColumn =settings.numVarsY;
    ydesignVar =   floor( (e-1)/numVarsinRow)+1;
    xdesignVar = mod(e-1,numVarsinRow)+1;
    
    % trying to find the nodes and then dofs of the macro problem
    shiftElementsX = (xdesignVar-1)*settings.numXElmPerDV;
    elementsInRow = settings.nelx;
    shiftNodesY = (ydesignVar-1)*settings.numYElmPerDV*elementsInRow;
    
%     count1 = 1;
%     count2 = 1;
%     isfirstrow = 1;
    
    
    nodelist=[];
    for i  =1: settings.numYElmPerDV
        startElement = shiftElementsX+shiftNodesY+(i-1)*elementsInRow;
        for j =1: settings.numXElmPerDV
            ee = startElement+j;
            
            nodes1=  IEN(ee,:);
            nodelist = [nodelist nodes1];
        end
    end
    
    nodelist = unique(sort(nodelist));
    
    % nodes1=  IEN(e,:);
    macroElementProps.elementNodes=nodelist;
%     [t1,t2]= size(nodelist);
    
    
    % multiloading cases. 
     [~, t2] = size(settings.loadingCase);        
    for loadcaseIndex = 1:t2
         utemp = U(loadcaseIndex,:);   
        
        displacementsX = utemp(nodelist*2-1);
        displacementsY = utemp(nodelist*2);

        macroElementProps.xDisplacements(loadcaseIndex,:)=displacementsX;
        macroElementProps.yDisplacements(loadcaseIndex,:)=displacementsY;
    end

    macroElementProps.elementNumber = e;
    % Save element to XY position map (needed for x and w vars retrival)
    %         outname = sprintf('./out%i/elementXYposition%i.csv',folderNum,macro_meso_iteration);
    %         elementXYposition=csvread(outname);
    %         results = elementXYposition(macroElementProps.elementNumber,:);
     macroElementProps.yPosition = ydesignVar;
     macroElementProps.xPosition = xdesignVar;
    
    % Get the density field
    outname = sprintf('./out%i/densityfield%i.csv',folderNum,macro_meso_iteration);
    x = csvread(outname);
    macroElementProps.density = x(macroElementProps.yPosition,macroElementProps.xPosition );
    
    % % Get the volume fraction field
    outname = sprintf('./out%i/volfractionfield%i.csv',folderNum,macro_meso_iteration);
    w = csvread(outname);
    macroElementProps.material1Fraction = w(macroElementProps.yPosition,macroElementProps.xPosition );
    
    % if(macro_meso_iteration>1)
    outname = sprintf('./out%i/Dgiven_%i_forElement_%i.csv',folderNum,macro_meso_iteration,e);
    D = csvread(outname);
    macroElementProps.D_given =D;
    % end
    
    
    
end