classdef DesignVars
    % Design varriables and support temp varriables are stored in this class.
    
    properties
        % --------------------------------
        % Design Var Arrays
        % --------------------------------
        x; % the "density" at each element
        w; % the volume fraction at each element
        Exx; % the Exx of the orthogonal material
        Eyy; % the Exx of the orthogonal material
        t; %theta, the rotation of the orthogonal material
        
        % Sub System Copies of design Vars and Lagrangian Multipliers
        lambdaExx=0; % lambda for the Exx
        lambdaEyy=0; % lambda for the Eyy
        lambdaTheta=0; % lambda for the Theta
        penaltyExx=0; % penality value used for the Exx Augmented Lagrangian Multiplier
        penaltyEyy=0; % penality value used for the  Eyy Augmented Lagrangian Multiplier
        penaltyTheta=0; %penality value used for the Theta Augmented Lagrangian Multiplier
        
        ExxSub=0; % Sub system copies of design var
        EyySub=0; % Sub system copies of design var
        thetaSub=0; % Sub system copies of design var
        
        
        
        % Optimization vars (volume fraction optimization)
        lambda1 = 0;
        mu1 = 1;
        c = 0; % objective.
        cCompliance = 0;
        cHeat = 0;
        storeOptimizationVar= []
        
        % --------------------------
        % Support vars
        % -------------------------------
        %         xold; %
        sensitivityElastic; % Sensitivity 1, elastic
        sensitivityElasticPart2; % Needed for E_yy optimization
        sensitivityHeat; % Sensitivity 2, heat
        currentVol1Fraction;
        currentVol2Fraction;
        targetAverageE =0;
        actualAverageE=0;
        
        %complianceSensitivity; %
        totalStress;
        dc; % Derivative of c (hence dc). C is the objective.
        dcSum; % derivitive of meso element sensitivies when designing a single meso structure for the whole macro structure.
        %g1elastic; % Derivative of c with respect to a material change for the elastic
        %g1heat; %  Derivative of cHeat with respect to a material change for the heat
        IEN; % element to node map. Save this matrix, so it does not need to be recalculated every time.
        XLocations; %=zeros(numNodesInRow,numNodesInColumn);
        YLocations; %=zeros(numNodesInRow,numNodesInColumn);
        globalPosition; %  = zeros(nn,2);
        NodeToXYArrayMap; % map of node numbers to their X,Y position in FEA arrays
        elementXYposition; % or othe position in the X ad W matrixes
        
        U_heatColumn; % temperature matrix, gives the temperature at each node (not at the element centers)
        U; % displacement matrix
        B; %  holding the displacement strain matrix
        
        maxF = 0;
        maxU = 0;
        
        
        % meso informatoin
        nelxTile;
        nelyTile
        xTile;
        wTile;
        IENTile;
        XLocationsTile;
        YLocationsTile;
        globalPositionTile;
        NodeToXYArrayMapTile;
        UTile;
        d11;
        d12;
        d22;
        d33;
        De11
        De12
        De22
        De33;
        
        mesoAddAdjcentCellDataObject;
        
        
        % When multiple elements are controlled by a single design var,
        % this array gives you the design var number.
        % Give it a X,Y position and it will tell you the the rows are the
        % element numbers.
        elementToDesignVarMap;
    end
    
    methods
        % Constructur method
        function obj = DesignVars(config)
            
            % Get the B matrix, the E,v,G do not matter and are not used
            % in the B calculation, so set them to 1.
            E = 1; v= 1; G = 1; strain = [];
            [~, ~, B_out] = elK_elastic(E,v, G,strain,[]);
            obj.B = B_out;
            
            macroDesignMode = 0;
            if(config.mode<100)
                macroDesignMode=1;
            end
            
            if ( macroDesignMode==1)
                % ------------------------------
                % Initialize Arrays used for Macro Optimization
                % ------------------------------
                obj.x(1:config.nely,1:config.nelx) =config.v1+config.v2; % artificial density of the elements
                obj.w(1:config.nely,1:config.nelx)  =config.v1; % actual volume fraction composition of each element
                %                   obj.w(1:config.nely,1:config.nelx)  =ones(config.nely,config.nelx);
                %                 obj.d(1:config.nely,1:config.nelx)  =ones(config.nely,config.nelx)*0.8; % orthotropic masterial distribution
                obj.Exx(1:config.nely,1:config.nelx)  =ones(config.nely,config.nelx);
                obj.Eyy(1:config.nely,1:config.nelx)  =ones(config.nely,config.nelx);
                obj.t(1:config.nely,1:config.nelx)  =ones(config.nely,config.nelx)*0; % rotation of the orthotropic material
                
                obj.sensitivityElastic(1:config.nely,1:config.nelx) = 0;
                obj.sensitivityElasticPart2 (1:config.nely,1:config.nelx) = 0;
                obj.sensitivityHeat(1:config.nely,1:config.nelx) = 0;
                if (config.doPlotStress == 1)
                    obj.totalStress(1:config.nely,1:config.nelx) = 0;
                end
                obj=obj.CalcIENmatrix(config);
                obj=obj.CalcNodeLocation(config);
                obj=obj.PreCalculateXYmapToNodeNumber(config);
                obj=obj.CalcElementXYposition(config);
                %     matProp=  matProp.ReadConstitutiveMatrixesFromFiles(config);
                
                
                % Start them out as 0. Set later if iteration 2 or higher
                obj.lambdaExx= obj.Exx*0; % lambda for the Exx
                obj.lambdaEyy= obj.Exx*0; % lambda for the Eyy
                obj.lambdaTheta=obj.Exx*0; % lambda for the Theta
                obj.penaltyExx= obj.Exx*0; % penality value used for the Exx Augmented Lagrangian Multiplier
                obj.penaltyEyy= obj.Exx*0; % penality value used for the  Eyy Augmented Lagrangian Multiplier
                obj.penaltyTheta= obj.Exx*0; %penality value used for the Theta Augmented Lagrangian Multiplier
                
                obj.ExxSub=obj.Exx*0; % Sub system copies of design var
                obj.EyySub=obj.Exx*0; % Sub system copies of design var
                obj.thetaSub=obj.Exx*0; % Sub system copies of design var
                
                
                
            end
            
            
            
        end
        
        % -----------------------------
        % Get the state of the macro optimization saved in csv files.
        %
        % Also, read Exx,Eyy,Theta, Lambdas and ExxSub, EyySub, ThetaSub
        % -----------------------------
        function obj = GetMacroStateVarsFromCSVFiles(obj, config)
            
            if(config.macro_meso_iteration>1)
                
                if (config.mode==60  )
                    % -------------------------
                    % Read from csv
                    % x (topology var)
                    % Exx System
                    % Eyy System
                    % t (Theta) System
                    %
                    % lambdaExx Array
                    % lambdaEyy Array
                    % lambdaTheta Array
                    %
                    % Exx SubSystem
                    % Eyy SubSystem
                    % t (Theta) SubSystem
                    %
                    % penality Exx
                    % penalty Eyy
                    % penalty theta
                    %
                    % See function SaveMacroProblemStateToCSV and GenerateRhoFunctionOfExxEyy
                    % -------------------------
                    folderNum = config.iterationNum;
                    oldIteration = config.macro_meso_iteration-1;
                    
                    outnameX = sprintf('./out%i/SIMPdensityfield%i.csv',folderNum,oldIteration);
                    outnameExx = sprintf('./out%i/ExxValues%i.csv',folderNum,oldIteration);
                    outnameEyy = sprintf('./out%i/EyyValues%i.csv',folderNum,oldIteration);
                    outnameTheta = sprintf('./out%i/ThetaValues%i.csv',folderNum,oldIteration);
                    outnamelambdaExx = sprintf('./out%i/lambdaExx%i.csv',folderNum,oldIteration);
                    outnamelambdaEyy = sprintf('./out%i/lambdaEyy%i.csv',folderNum,oldIteration);
                    outnamelambdaTheta = sprintf('./out%i/lambdaTheta%i.csv',folderNum,oldIteration);
                    outnamepenaltyExx = sprintf('./out%i/penaltyExx%i.csv',folderNum,oldIteration);
                    outnamepenaltyEyy = sprintf('./out%i/penaltyEyy%i.csv',folderNum,oldIteration);
                    outnamepenaltyTheta = sprintf('./out%i/penaltyTheta%i.csv',folderNum,oldIteration);
                    outnameExxSubSysValues = sprintf('./out%i/ExxSubSysValues%i.csv',folderNum,oldIteration);
                    outnameEyySubSysValues = sprintf('./out%i/EyySubSysValues%i.csv',folderNum,oldIteration);
                    outnameThetaSubSysValues = sprintf('./out%i/ThetaSubSysValues%i.csv',folderNum,oldIteration);
                    
                    
                    obj.x = csvread(outnameX); % the "density" at each element
                    obj.Exx = csvread(outnameExx); % the Exx of the orthogonal material
                    obj.Eyy = csvread(outnameEyy); % the Eyy of the orthogonal material
                    obj.t = csvread(outnameTheta); %theta, the rotation of the orthogonal material
                    
                    % Sub System Copies of design Vars and Lagrangian Multipliers
                    % starting on iteration 3, get the saved values.
                    % these values are initilized on iteration 2
                    if(config.macro_meso_iteration>2)
                        obj.lambdaExx= csvread(outnamelambdaExx); % lambda for the Exx
                        obj.lambdaEyy= csvread(outnamelambdaEyy); % lambda for the Eyy
                        obj.lambdaTheta= csvread(outnamelambdaTheta); % lambda for the Theta
                        obj.penaltyExx= csvread(outnamepenaltyExx); % penality value used for the Exx Augmented Lagrangian Multiplier
                        obj.penaltyEyy= csvread(outnamepenaltyEyy); % penality value used for the  Eyy Augmented Lagrangian Multiplier
                        obj.penaltyTheta= csvread(outnamepenaltyTheta); %penality value used for the Theta Augmented Lagrangian Multiplier
                    elseif(config.macro_meso_iteration==2)
                        % set the lagrangians to zero for now. They will be
                        % updated later after the penalty is calculated. 
                        obj.lambdaExx=zeros(config.nely,config.nelx);
                        obj.lambdaEyy= obj.lambdaExx;
                        obj.lambdaTheta= obj.lambdaExx;                     
                    end
                    
                    obj.ExxSub=csvread(outnameExxSubSysValues); % Sub system copies of design var
                    obj.EyySub=csvread(outnameEyySubSysValues); % Sub system copies of design var
                    obj.thetaSub=csvread(outnameThetaSubSysValues); % Sub system copies of design var
                    
                    
                elseif(config.mode==50)
                    % -------------------------
                    % (old) if Volume Fraction optimization
                    % -------------------------
                    folderNum = config.iterationNum;
                    previousIterationNum = config.macro_meso_iteration-1;
                    
                    % get the topology densities.
                    outname = sprintf('./out%i/densityfield%i.csv',folderNum,previousIterationNum);
                    obj.x = csvread(outname);
                    
                    % get the volume fraction optimization vars
                    outname = sprintf('./out%i/volfractionfield%i.csv',folderNum,previousIterationNum);
                    obj.w=csvread(outname);
                    
                    % get the lambda value
                    outname = sprintf('./out%i/lambda%i.csv',folderNum,previousIterationNum);
                    obj.lambda1=csvread(outname);
                else
                    message = 'First iteration, no macro state vars to read';
                end
            end
        end
        
        
        
        % -----------------------------
        % UpdatePenaltyAndLagrangianValues
        % -----------------------------
        function obj = UpdatePenaltyAndLagrangianValues(obj, config,matProp)
            
             ne = config.nelx*config.nely;
                  diffExx = obj.Exx-obj.ExxSub;
                diffEyy = obj.Eyy - obj.EyySub;
                diffTheta = obj.t-obj.thetaSub;
             
            % if the first time, calculate the initial penalty values. 
            if(config.macro_meso_iteration==2)
                % Get displacement fields
                oldIteration = config.macro_meso_iteration-1;
                folderNum = config.iterationNum;
                outname = sprintf('./out%i/displacement%i.csv',folderNum,oldIteration);
                UpreviousIteration =  csvread(outname);
                
           
                
                % loop over the elements.
               
                for e = 1:ne
                    
                    [xPos,yPos]= obj.GivenNodeNumberGetXY(e);
                    
                    nodes1 = obj.IEN(e,:);
                    xNodes = nodes1*2-1;
                    yNodes = nodes1*2;
                    
                    % I cannot use the union, or else the order get messed up. The order
                    % is important. Same in the actual topology code when you are
                    % calculating the objectiv
                    NodeNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
                    
                    
                    KE = matProp.getKMatrixTopExxYyyRotVars(config,obj.x(yPos,xPos),obj.Exx(yPos,xPos), obj.Eyy(yPos,xPos),obj.t(yPos,xPos),obj.w(yPos,xPos));
                    
                    
                    % get the displacements and calculate the strain energy
                    strainEnergy = 0;
                    [~, t2] = size(config.loadingCase);
                    for loadcaseIndex = 1:t2
                        u_local = transpose(UpreviousIteration(loadcaseIndex,NodeNumbers));
                        strainEnergy_temp=transpose(u_local)*KE*u_local;
                        strainEnergy=strainEnergy+strainEnergy_temp;
                      
                        
                    end
                   
                    % ---------------------------------------------
                    % Calculate the penalty values!!
                    % ---------------------------------------------
                    omegaLocal = config.Omega;
                    obj.penaltyExx(yPos,xPos)=2*omegaLocal*strainEnergy/ abs( diffExx(yPos,xPos));
                    obj.penaltyEyy(yPos,xPos)=2*omegaLocal*strainEnergy/ abs (diffEyy(yPos,xPos));
                    obj.penaltyTheta(yPos,xPos)=2*omegaLocal*strainEnergy/ abs( diffTheta(yPos,xPos));
                    
                end
            end
            
            
            
            % ---------------------------------------------
            % Calculate the penalty values!! if on iteration 3 or
            % greater
            % ---------------------------------------------
            if(config.macro_meso_iteration>2)
                
                updateMultiplier = 3;
                for e = 1:ne                    
                    [xPos,yPos]= obj.GivenNodeNumberGetXY(e);                    
                     obj.penaltyExx(yPos,xPos)=   obj.penaltyExx(yPos,xPos)*updateMultiplier;
                    obj.penaltyEyy(yPos,xPos)= obj.penaltyEyy(yPos,xPos)*updateMultiplier;
                    obj.penaltyTheta(yPos,xPos)=  obj.penaltyTheta(yPos,xPos)*updateMultiplier;                    
                end
                
            end
            
            
             % ---------------------------------------------
             % Calculate the lagrangian values.
             % ---------------------------------------------
             
             for e = 1:ne
                 [xPos,yPos]= obj.GivenNodeNumberGetXY(e);
                 obj.lambdaExx(yPos,xPos)= obj.lambdaExx(yPos,xPos)+obj.penaltyExx(yPos,xPos) *diffExx(yPos,xPos);
                 obj.lambdaEyy(yPos,xPos)=  obj.lambdaEyy(yPos,xPos)+obj.penaltyEyy(yPos,xPos) *diffEyy(yPos,xPos);
                 obj.lambdaTheta(yPos,xPos)=   obj.lambdaTheta(yPos,xPos)+obj.penaltyTheta(yPos,xPos) *diffTheta(yPos,xPos);
             end
             
             if 1==2
                 p = plotResults;
                 figure
                 subplot(2,3,1)
                 p.PlotArrayGeneric( obj.penaltyExx, ' obj.penaltyExx')
                 subplot(2,3,2)
                 p.PlotArrayGeneric( obj.penaltyEyy, ' obj.penaltyEyy')
                 subplot(2,3,3)
                 p.PlotArrayGeneric( obj.penaltyTheta, ' obj.penaltyTheta')
                 
                 subplot(2,3,4)
                 p.PlotArrayGeneric( obj.lambdaExx, ' obj.lambdaExx')
                 subplot(2,3,5)
                 p.PlotArrayGeneric( obj.lambdaEyy, ' obj.lambdaEyy')
                 subplot(2,3,6)
                 p.PlotArrayGeneric( obj.lambdaTheta, ' obj.lambdaTheta')
             end
        end
        
        
        
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
        % dof
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
        function obj =  CalcIENmatrix(obj,config)
            
            count = 1;
            elementsInRow = config.nelx+1; % think this is actually "nodes in a row. "
            nn = (config.nelx+1)*(config.nely+1); % number of nodes
            obj.IEN = zeros(nn,4);
            % Each row, so nely # of row
            for i = 1:config.nely
                rowMultiplier = i-1;
                % Each column, so nelx # of row
                for j= 1:config.nelx
                    obj.IEN(count,:)=[rowMultiplier*elementsInRow+j, ...
                        rowMultiplier*elementsInRow+j+1, ...
                        (rowMultiplier +1)*elementsInRow+j+1,...
                        (rowMultiplier +1)*elementsInRow+j];
                    count = count+1;
                end
            end
        end
        
        
        % Calcualte the node locations and put the information
        % into an array. Needed for the FEA
        % Calculate it here, so it only need to be calculated once.
        %
        % Think these are actually node locations
        function obj = CalcNodeLocation(obj,config)
            nn = (config.nelx+1)*(config.nely+1); % number of nodes
            
            numNodesInRow = config.nelx + 1;
            numNodesInColumn = config.nely + 1;
            %             obj.XLocations=zeros(numNodesInRow,numNodesInColumn);
            %             obj.YLocations=zeros(numNodesInRow,numNodesInColumn);
            
            obj.globalPosition = zeros(nn,2);
            count = 1;
            for i = 1:numNodesInColumn  % y
                for j= 1:numNodesInRow % x
                    obj.globalPosition(count,:) = [j-1 i-1];
                    count = count+1;
                    %                     obj.XLocations(j,i) = j-1;
                    %                     obj.YLocations(j,i) = i-1;
                end
            end
        end
        
        % Pre Calculate a map of the XY array coordinates for each node
        % number.
        function obj= PreCalculateXYmapToNodeNumber(obj ,config)
            nn = (config.nelx+1)*(config.nely+1); % number of nodes
            obj.NodeToXYArrayMap = zeros(nn,2);
            count = 1;
            for i = 1:config.nely
                for j= 1:config.nelx
                    obj.NodeToXYArrayMap(count,:) = [i,j];
                    count = count+1;
                end
            end
        end
        
        % -----------------------------
        %
        % Get the design var position when given an element X,Y position.
        %
        % Only applicable when multiple elements per design var is true.
        %
        % -----------------------------
        %         function [x,y] = GetDesignVarPositionGivenXYElement(obj,config,xelm,yelm)
        %             x = ceil(xelm/config.numXElmPerDV);
        %             y = ceil(yelm/config.numYElmPerDV);
        %         end
        
        
        
        % -----------------------------
        % Calcualte the Center of each element and put the information
        % into an array. Needed for the FEA
        % Calculate it here, so it only need to be calculated once.
        % -----------------------------
        function obj = CalcNodeLocationMeso(obj,config)
            nn = (config.nelx)*(config.nely); % number of nodes
            
            numNodesInRow = config.nelx ;
            numNodesInColumn = config.nely ;
            obj.XLocations=zeros(numNodesInRow,numNodesInColumn);
            obj.YLocations=zeros(numNodesInRow,numNodesInColumn);
            
            obj.globalPosition = zeros(nn,2);
            count = 1;
            for i = 1:(numNodesInColumn)  % y
                for j= 1:(numNodesInRow) % x
                    obj.globalPosition(count,:) = [j-1 i-1];
                    count = count+1;
                    obj.XLocations(j,i) = j-1;
                    obj.YLocations(j,i) = i-1;
                end
            end
        end
        
        % -----------------------------
        % -----------------------------
        function obj = CalcNodeLocationMeso_Tile(obj,config)
            nn = (obj.nelxTile)*(obj.nelyTile); % number of nodes same as number of elements since it wraps
            %             nn = (config.nelx)*(config.nely); % number of nodes
            
            numNodesInRow = obj.nelxTile ;
            numNodesInColumn = obj.nelyTile;
            obj.XLocations=zeros(numNodesInRow,numNodesInColumn);
            obj.YLocations=zeros(numNodesInRow,numNodesInColumn);
            
            obj.globalPosition = zeros(nn,2);
            count = 1;
            for i = 1:(numNodesInColumn)  % y
                for j= 1:(numNodesInRow) % x
                    obj.globalPositionTile(count,:) = [j-1 i-1];
                    count = count+1;
                    obj.XLocationsTile(j,i) = j-1;
                    obj.YLocationsTile(j,i) = i-1;
                end
            end
        end
        
        
        
        % ---------------------------
        %
        %  Calcualte the element node map (IEN) for the homogenization meso
        %  structure case. BAsically the map needs to have the displacement
        %  field be periodic, so ti must loop back on itself on the edges.
        %
        % ---------------------------
        function obj =  CalcElementNodeMapmatrixWithPeriodicXandY(obj,config)
            
            count = 1;
            elementsInRow = config.nelx;
            nn = (config.nelx)*(config.nely); % number of nodes same as number of elements since it wraps
            obj.IEN = zeros(nn,4);
            % Each row, so nely # of row
            for i = 1:config.nely
                rowMultiplier = i-1;
                % Each column, so nelx # of row
                for j= 1:config.nelx
                    
                    % normal case
                    if(j ~= config.nelx && i ~= config.nely )
                        obj.IEN(count,:)=[rowMultiplier*elementsInRow+j, ...
                            rowMultiplier*elementsInRow+j+1, ...
                            (rowMultiplier +1)*elementsInRow+j+1,...
                            (rowMultiplier +1)*elementsInRow+j];
                    elseif(j == config.nelx && i ~= config.nely )
                        % On the right side of the mesh case
                        
                        % loop back arround case case
                        obj.IEN(count,:)=[rowMultiplier*elementsInRow+j, ...
                            rowMultiplier*elementsInRow+1, ... % note the difference
                            (rowMultiplier +1)*elementsInRow+1,... % note the difference
                            (rowMultiplier +1)*elementsInRow+j];
                    elseif(j ~= config.nelx && i == config.nely )
                        % On the top side of the mesh case
                        
                        obj.IEN(count,:)=[rowMultiplier*elementsInRow+j, ...
                            rowMultiplier*elementsInRow+j+1, ...
                            0+j+1,... % note the difference
                            0+j]; % note the difference
                        
                    elseif(j == config.nelx && i == config.nely )
                        % On the top right  side of the mesh case
                        obj.IEN(count,:)=[rowMultiplier*elementsInRow+j, ...
                            rowMultiplier*elementsInRow+1, ...
                            1,... % back to the node 1
                            j];
                        
                    end
                    count = count+1;
                end
            end
        end
        
        
        
        
        % --------------------------------------------
        %
        %       ELEMENT TO NODE MAP FOR TILED PERIODIC MESO
        %
        % --------------------------------------------
        function obj =  CalcElementNodeMapmatrixWithPeriodicXandY_Tile(obj,config)
            
            count = 1;
            elementsInRow=obj.nelxTile;% elementsInRow = config.nelx;
            nn = (obj.nelxTile)*(obj.nelyTile); % number of nodes same as number of elements since it wraps
            obj.IENTile = zeros(nn,4);
            % Each row, so nely # of row
            for i = 1:obj.nelyTile
                rowMultiplier = i-1;
                % Each column, so nelx # of row
                for j= 1:obj.nelxTile
                    
                    % normal case
                    if(j ~= obj.nelxTile && i ~= obj.nelyTile )
                        obj.IENTile(count,:)=[rowMultiplier*elementsInRow+j, ...
                            rowMultiplier*elementsInRow+j+1, ...
                            (rowMultiplier +1)*elementsInRow+j+1,...
                            (rowMultiplier +1)*elementsInRow+j];
                    elseif(j == obj.nelxTile && i ~= obj.nelyTile )
                        % On the right side of the mesh case
                        
                        % loop back arround case case
                        obj.IENTile(count,:)=[rowMultiplier*elementsInRow+j, ...
                            rowMultiplier*elementsInRow+1, ... % note the difference
                            (rowMultiplier +1)*elementsInRow+1,... % note the difference
                            (rowMultiplier +1)*elementsInRow+j];
                    elseif(j ~= obj.nelxTile && i == obj.nelyTile )
                        % On the top side of the mesh case
                        
                        obj.IENTile(count,:)=[rowMultiplier*elementsInRow+j, ...
                            rowMultiplier*elementsInRow+j+1, ...
                            0+j+1,... % note the difference
                            0+j]; % note the difference
                        
                    elseif(j == obj.nelxTile && i == obj.nelyTile )
                        % On the top right  side of the mesh case
                        obj.IENTile(count,:)=[rowMultiplier*elementsInRow+j, ...
                            rowMultiplier*elementsInRow+1, ...
                            1,... % back to the node 1
                            j];
                        
                    end
                    count = count+1;
                end
            end
        end
        
        % Given an X, find the node number
        function number = GetNodeNumberGivenXY(config, x,y)
            numNodesInRow = config.nelx+1;
            % numNodesInColumn = obj.nely+1;
            rowMultiplier = y-1;
            number = rowMultiplier*numNodesInRow+x;
        end
        
        % Given a node number, find the X, Y position (not the physical
        % position, but the matrix location position)
        function [x , y ]= GivenNodeNumberGetXY(obj, nodeNum)
            [result] =  obj.NodeToXYArrayMap(nodeNum,:);
            y = result(1);
            x = result(2);
        end
        
        % Calculates the volume
        function [obj] = CalculateVolumeFractions(obj, config,matProp)
            
            
            %             Xtemp = obj.x;
            %             Xtemp(obj.x>config.voidMaterialDensityCutOff)=1;
            %             Xtemp(obj.x<=config.voidMaterialDensityCutOff)=0;
            %             ne = config.nelx*config.nely;
            %                ne = config.nelx*config.nely;
            %                 averageElasticLocal = (sum(sum(obj.Exx.*obj.x))+sum(sum(obj.Eyy.*obj.x)))/ne;
            %                 volume1 = sum(sum(obj.x))/ne;
            %                 volume2 = averageElasticLocal/matProp.E_material1;
            %                   volume1 = (sum(sum(obj.Exx.*obj.x))+sum(sum(obj.Eyy.*obj.x)))/ne;
            
            ne = config.nelx*config.nely;
            %             neSolid = config.nelx*config.nely*(config.v1+config.v2);
            totalMaterial = sum(sum(obj.x));
            obj.targetAverageE=(config.v1*matProp.E_material1+config.v2*matProp.E_material2)/(config.v1+config.v2);
            
            if(config.useExxEyy==1)
                %                    avg= 0.5*(obj.Exx+obj.Eyy);
                totalExx =obj.x.*obj.Exx;
                totalEyy = obj.x.* obj.Eyy;
                avgE = (totalExx+totalEyy)/2;
                if(config.testingVerGradMaterail ==1)
                    minE = matProp.E_material2;
                else
                    minE = matProp.E_material2/2;
                end
                obj.w = avgE/(matProp.E_material1-minE);
                %                   for i = 1:config.nelx
                %                       for j = 1:config.nely
                %                            obj.w(j,i) =max(obj.Exx(j,i),obj.Eyy(j,i))/(matProp.E_material1-minE);
                %                       end
                %                   end
                %                   obj.w =max(obj.Exx,obj.Eyy)/matProp.E_material1;
                
                
                obj.actualAverageE= sum(sum(avgE))/totalMaterial;
                
                % volume fraction is over the domain, while target E
                obj.  currentVol1Fraction =sum(sum( obj.x.*obj.w))/ne;
                obj.   currentVol2Fraction =sum(sum( obj.x.*(1-obj.w)))/ne;
            else
                
                
                
                totalMat1 =sum(sum( obj.x.*obj.w*matProp.E_material1));
                totalMat2 =sum(sum( obj.x.*(1-obj.w)*matProp.E_material2));
                % obj.actualAverageE= obj.currentVol1Fraction*matProp.E_material1+  obj. currentVol2Fraction*matProp.E_material2;
                obj.actualAverageE= (totalMat1+totalMat2)/totalMaterial;
                obj.  currentVol1Fraction =sum(sum( obj.x.*obj.w))/ne;
                obj.   currentVol2Fraction =sum(sum( obj.x.*(1-obj.w)))/ne;
            end
            
            
        end
        
        
        %---------------------------
        % Combine the Exx and Eyy vars into the vol fraction (w) var.
        % This is just used so that we can plot more easily.
        %---------------------------
        %         function [obj] = CalcVolFractionUsingExxEyy(obj, config, matProp)
        %             avg= 0.5*(obj.Exx+obj.Eyy);
        %
        %             minE = matProp.E_material1/2;
        %             temp = avg-minE;
        %
        %             w = temp/(matProp.E_material1-minE);
        %
        %              obj.w = w;
        %         end
        
        %%
        
        % --------------------------------------------
        %
        %      Calculate the element XY position.
        %
        % --------------------------------------------
        function obj = CalcElementXYposition(obj,config)
            nelm = config.nelx*config.nely;
            
            obj.elementXYposition =   zeros(nelm,2);
            
            e_count = 1;
            for ely = 1:config.nely
                %                 rowMultiplier = ely-1;
                for elx = 1:config.nelx
                    obj.elementXYposition(e_count,:)=[ely,elx];
                    e_count = e_count+1;
                end
            end
        end
        
        
        % ----------------------------------
        %
        %  RUN THE FEA ANAYLSIS FOR EACH LOADING CASE
        %
        % --------------------------------
        function obj = RunFEAs(obj, config, matProp, loop)
            [~, t2] = size(config.loadingCase);
            for loadcaseIndex = 1:t2
                loadcase = config.loadingCase(loadcaseIndex);
                
                % FE-ANALYSIS
                if (config.w1 ~= 1)
                    u_heat_loadcase  =temperatureFEA_V3(obj, config, matProp,loop,loadcase);
                    obj.U_heatColumn(loadcaseIndex,:)=u_heat_loadcase;
                end
                [UloadCase, obj.maxF, obj.maxU]=FE_elasticV2(obj, config, matProp,loadcase);
                obj.U(loadcaseIndex,:)=UloadCase;
            end
            
            %             obj.c = 0.; % c is the objective. Total strain energy
            %             obj.cCompliance = 0;
            %             obj.cHeat = 0;
        end
        
        % -----------------------------------------
        %
        % Calculate Topology Sensitivity
        %
        % -- FEA data sent to this function.
        %
        % -----------------------------------------
        function obj = CalculateTopologySensitivity(obj, config, matProp, loop)
            elementsInRow = config.nelx+1;
            
            obj.sensitivityHeat(:,:)=0;
            obj.sensitivityElastic(:,:)=0;
            obj.cCompliance=0;
            
            % allow multiple loading cases.
            [~, t2] = size(config.loadingCase);
            
            %             for loadcaseIndex = 1:t2
            % UloadCase= obj.U(loadcaseIndex,:);
            
            % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
            count =1;
            for ely = 1:config.nely
                rowMultiplier = ely-1;
                for elx = 1:config.nelx
                    nodes1=[rowMultiplier*elementsInRow+elx;
                        rowMultiplier*elementsInRow+elx+1;
                        (rowMultiplier +1)*elementsInRow+elx+1;
                        (rowMultiplier +1)*elementsInRow+elx];
                    
                    xNodes = nodes1*2-1;
                    yNodes = nodes1*2;
                    NodeNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
                    
                    % if heat objective!!!
                    if (config.w1 ~= 1)
                        U_heat = obj.U_heatColumn(nodes1,:);
                        %                             averageElementTemp = mean2(U_heat); % calculate the average temperature of the 4 nodes
                        KEHeat = matProp.effectiveHeatKEmatrix(  obj.w(ely,elx), config);
                        obj.cHeat =   obj.cHeat           + obj.x(ely,elx)^config.penal*U_heat'*KEHeat*U_heat;
                        
                        % calculate the minim temp sensitivity
                        obj.sensitivityHeat(ely,elx) = -config.penal*obj.x(ely,elx)^(config.penal-1)*U_heat'*KEHeat*U_heat + obj.sensitivityHeat(ely,elx);
                        % obj.g1heat(ely,elx) = obj.x(ely,elx)^(config.penal)*U_heat'*matProp.dKheat*U_heat + obj.g1heat(ely,elx) ;
                    end
                    
                    % Calculate generic KE with 1 as the SIMP density
                    KEgeneric = matProp.getKMatrixTopExxYyyRotVars(config,1,obj.Exx(ely,elx), obj.Eyy(ely,elx),obj.t(ely,elx),obj.w(ely,elx));
                    % Calculate the elastic sensitivity!
                    KEsensitive = KEgeneric*config.penal*obj.x(ely,elx)^(config.penal-1);
                    % Calculate the actual KE
                    KE = KEgeneric*obj.x(ely,elx)^(config.penal);
                    
                    % allow multiple loading cases.
                    Urows = obj.U(:,NodeNumbers);
                    for tt = 1:t2
                        Ue = Urows(tt,:)';
                        obj.sensitivityElastic(ely,elx) =-Ue'*KEsensitive*Ue+ obj.sensitivityElastic(ely,elx);
                        % Calculate the total elastic compliance
                        obj.cCompliance = obj.cCompliance + Ue'*KE*Ue;
                    end
                    
                    count=count+1;
                    
                end % end, loop over x
            end % end, for loop over nely
            %             end % end for loop over load cases.
            
            obj.c=obj.cCompliance*config.w1+obj.cHeat*config.w2;
        end % End Function, CalculateSenstivities
        
        % -----------------------------------------
        %
        % CalculateObjectiveValue
        %
        %  % -- FEA data sent to this function.
        %
        % -----------------------------------------
        function obj = CalculateObjectiveValue(obj, config, matProp, loop)
            elementsInRow = config.nelx+1;
            
            obj.cCompliance=0;
            
            % allow multiple loading cases.
            [~, t2] = size(config.loadingCase);
            
            %             for loadcaseIndex = 1:t2
            % UloadCase= obj.U(loadcaseIndex,:);
            
            % OBJECTIVE FUNCTION
            count =1;
            for ely = 1:config.nely
                rowMultiplier = ely-1;
                for elx = 1:config.nelx
                    nodes1=[rowMultiplier*elementsInRow+elx;
                        rowMultiplier*elementsInRow+elx+1;
                        (rowMultiplier +1)*elementsInRow+elx+1;
                        (rowMultiplier +1)*elementsInRow+elx];
                    
                    xNodes = nodes1*2-1;
                    yNodes = nodes1*2;
                    NodeNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
                    
                    % if heat objective!!!
                    if (config.w1 ~= 1)
                        U_heat = obj.U_heatColumn(nodes1,:);
                        %                             averageElementTemp = mean2(U_heat); % calculate the average temperature of the 4 nodes
                        KEHeat = matProp.effectiveHeatKEmatrix(  obj.w(ely,elx), config);
                        obj.cHeat =   obj.cHeat           + obj.x(ely,elx)^config.penal*U_heat'*KEHeat*U_heat;
                        
                        % calculate the minim temp sensitivity
                        obj.sensitivityHeat(ely,elx) = -config.penal*obj.x(ely,elx)^(config.penal-1)*U_heat'*KEHeat*U_heat + obj.sensitivityHeat(ely,elx);
                        % obj.g1heat(ely,elx) = obj.x(ely,elx)^(config.penal)*U_heat'*matProp.dKheat*U_heat + obj.g1heat(ely,elx) ;
                    end
                    
                    % Calculate generic KE with 1 as the SIMP density
                    KEgeneric = matProp.getKMatrixTopExxYyyRotVars(config,1,obj.Exx(ely,elx), obj.Eyy(ely,elx),obj.t(ely,elx),obj.w(ely,elx));
                    % Calculate the elastic sensitivity!
                    %                     KEsensitive = KEgeneric*config.penal*obj.x(ely,elx)^(config.penal-1);
                    % Calculate the actual KE
                    KE = KEgeneric*obj.x(ely,elx)^(config.penal);
                    
                    % allow multiple loading cases.
                    Urows = obj.U(:,NodeNumbers);
                    for tt = 1:t2
                        Ue = Urows(tt,:)';
                        %                         obj.sensitivityElastic(ely,elx) =-Ue'*KEsensitive*Ue+ obj.sensitivityElastic(ely,elx);
                        % Calculate the total elastic compliance
                        obj.cCompliance = obj.cCompliance + Ue'*KE*Ue;
                    end
                    
                    count=count+1;
                    
                end % end, loop over x
            end % end, for loop over nely
            %             end % end for loop over load cases.
            
            obj.c=obj.cCompliance*config.w1+obj.cHeat*config.w2;
        end
        
        
        % -----------------------------------------
        %
        % CalculateTMaterialGradientSensitivity
        %
        %  % -- FEA data sent to this function.
        %
        % -----------------------------------------
        function obj = CalculateMaterialGradientSensitivity(obj, config, matProp, loop)
            elementsInRow = config.nelx+1;
            
            % allow multiple loading cases.
            [~, t2] = size(config.loadingCase);
            
            obj.sensitivityHeat(:,:)=0;
            obj.sensitivityElastic(:,:)=0;
            
            
            %             for loadcaseIndex = 1:t2
            %                 UloadCase=obj.U(loadcaseIndex,:);
            %                 count =1;
            
            for ely = 1:config.nely
                rowMultiplier = ely-1;
                for elx = 1:config.nelx
                    nodes1=[rowMultiplier*elementsInRow+elx;
                        rowMultiplier*elementsInRow+elx+1;
                        (rowMultiplier +1)*elementsInRow+elx+1;
                        (rowMultiplier +1)*elementsInRow+elx];
                    
                    xNodes = nodes1*2-1;
                    yNodes = nodes1*2;
                    NodeNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
                    
                    %                         Ue = UloadCase(NodeNumbers)';
                    
                    % Calculate the elastic sensitivity!
                    KEsensitive = matProp.getKMatrixGradientMaterialSensitivity(config,obj.x(ely,elx),obj.Exx(ely,elx), obj.Eyy(ely,elx),obj.t(ely,elx));
                    
                    % allow multiple loading cases.
                    Urows = obj.U(:,NodeNumbers);
                    for tt = 1:t2
                        Ue = Urows(tt,:)';
                        obj.sensitivityElastic(ely,elx) =Ue'*KEsensitive*Ue+ obj.sensitivityElastic(ely,elx);
                        % Calculate the total elastic compliance
                        %                         obj.cCompliance = obj.cCompliance + Ue'*KE*Ue;
                    end
                    
                    
                    %                         obj.sensitivityElastic(ely,elx) =Ue'*KEsensitive*Ue+ obj.sensitivityElastic(ely,elx);
                    
                    
                    % if heat objective!!!
                    if (config.w1 ~= 1)
                        U_heat = obj.U_heatColumn(nodes1,:);
                        %                             averageElementTemp = mean2(U_heat); % calculate the average temperature of the 4 nodes
                        KEHeat = matProp.effectiveHeatKEmatrix(  obj.w(ely,elx), config);
                        obj.cHeat =   obj.cHeat  + obj.x(ely,elx)^config.penal*U_heat'*KEHeat*U_heat;
                        
                        % calculate the minim temp sensitivity
                        %                             obj.temp2(ely,elx) = -config.penal*obj.x(ely,elx)^(config.penal-1)*U_heat'*KEHeat*U_heat + obj.temp2(ely,elx);
                        obj.sensitivityHeat(ely,elx) = obj.x(ely,elx)^(config.penal)*U_heat'*matProp.dKheat*U_heat + obj.sensitivityHeat(ely,elx) ;
                    end
                    
                    %                         count=count+1;
                end % end, loop over x
            end % end, for loop over y
            
            %             end % end for loop over load cases.
        end % End Function, CalculateMaterialGradientSensitivity
        
        
        
        
        % -----------------------------------------
        %
        % CalculateExxEyySensitivity
        %
        %  % -- FEA data sent to this function.
        %
        % -----------------------------------------
        function obj = CalculateExxEyySensitivity(obj, config, matProp, loop)
            elementsInRow = config.nelx+1;
            obj.sensitivityHeat(:,:)=0;
            obj.sensitivityElastic(:,:)=0;
            obj.sensitivityElasticPart2(:,:) =0;
            
            % allow multiple loading cases.
            [~, t2] = size(config.loadingCase);
            
            %             for loadcaseIndex = 1:t2
            %                 UloadCase= obj.U(loadcaseIndex,:);
            
            %  SENSITIVITY ANALYSIS
            for y = 1:config.nely
                rowMultiplier = y-1;
                for xx = 1:config.nelx
                    nodes1=[rowMultiplier*elementsInRow+xx;
                        rowMultiplier*elementsInRow+xx+1;
                        (rowMultiplier +1)*elementsInRow+xx+1;
                        (rowMultiplier +1)*elementsInRow+xx];
                    
                    xNodes = nodes1*2-1;
                    yNodes = nodes1*2;
                    NodeNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
                    
                    URows = obj.U(:,NodeNumbers);
                    
                    % Get the sensitivity K matrix
                    % Set Exx = 1, Eyy = 0 to get sensitivity
                    KExx = matProp.getKMatrixTopExxYyyRotVars(config,obj.x(y,xx),1, 0,obj.t(y,xx),[]);
                    
                    % Set Eyy= 1,  Eyy = 0,to get sensitivity
                    KEyy = matProp.getKMatrixTopExxYyyRotVars(config,obj.x(y,xx),0, 1,obj.t(y,xx),[]);
                    
                    
                    % allow multiple loading cases.
                    tempSensi = 0;
                    tempSensiPart2=0;
                    for i = 1:t2
                        Ucase = URows(i,:)';
                        tempSensi= tempSensi+Ucase'*KExx*Ucase;
                        tempSensiPart2 = tempSensiPart2++Ucase'*KEyy*Ucase;
                    end
                    
                    obj.sensitivityElastic(y,xx) = tempSensi + obj.sensitivityElastic(y,xx);
                    obj.sensitivityElasticPart2(y,xx) = tempSensiPart2 + obj.sensitivityElasticPart2(y,xx);
                    
                end % end, loop over x
            end % end, for loop over nely
            %             end % end for loop over load cases.
            
        end % End Function, CalculateOthogonalDistributionSensitivity
        
        % ------------------------------------------------------------
        %
        % CalculateSensitiviesMesoStructure no periodic
        %
        % ------------------------------------------------------------
        function obj = CalculateSensitiviesMesoStructureNoPeriodic(obj, config, matProp, loop,macroElemProps, U)
            
            doplot = 0;
            doplotfinal = 0;
            if(doplot ==1 || doplotfinal==1)
                
                p = plotResults;
                figure(1);
            end
            
            % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
            obj.c = 0.; % c is the objective. Total strain energy
            obj.cCompliance = 0;
            obj.cHeat = 0;
            
            [~, t2] = size(config.loadingCase);
            % allow multiple loading cases.
            
            
            
            for loadcaseIndex = 1:t2
                Ucase = U(:,loadcaseIndex);
                loadcase = config.loadingCase(loadcaseIndex);
                ne = config.nelx*config.nely; % number of elements
                for e = 1:ne
                    
                    % loop over local node numbers to get their node global node numbers
                    nodes1 = obj.IEN(e,:);
                    [elx,ely]= obj.GivenNodeNumberGetXY(e);
                    
                    xNodes = nodes1*2-1;
                    yNodes = nodes1*2;
                    dofNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
                    
                    Ue = Ucase(dofNumbers);
                    
                    % U_heat = obj.U_heatColumn(nodes1,:);
                    %averageElementTemp = mean2(U_heat); % calculate the average temperature of the 4 nodes
                    
                    % Get the element K matrix for this partiular element
                    KE = matProp.effectiveElasticKEmatrix(  obj.w(ely,elx),config,[]);
                    
                    % KEHeat = matProp.effectiveHeatKEmatrix(  obj.w(ely,elx), config);
                    % Dmaterial = matProp.calculateEffectiveConstitutiveEquation( obj.w(ely,elx), config);
                    %                 config.nelx
                    % Find the elastic strain
                    elasticStrain = obj.B*Ue;
                    
                    % term1 = transpose(Ue)*KE*Ue*obj.x(ely,elx)^(config.penal-1)*config.penal;
                    
                    term1 = transpose(Ue)*KE*Ue*obj.x(ely,elx)^(config.penal-1)*config.penal;
                    
                    %                      term1_method2 = (eye(3)-elasticStrain);
                    %  term2 = 0;
                    % term3= 0;
                    
                    % Sum the elastic compliance terms.
                    % total = (term1 + term2 + term3);
                    obj.temp1(ely,elx) = term1+obj.temp1(ely,elx);
                    
                    if(doplot ==1)
                        %                     if(mod(e,10) ==0)
                        
                        p.PlotArrayGeneric(obj.temp1, 'plotting sensitivities while running. ')
                        drawnow
                        %                     end
                    end
                    % calculate the minim temp sensitivity
                    % obj.temp2(ely,elx) = -config.penal*obj.x(ely,elx)^(config.penal-1)*U_heat'*KEHeat*U_heat;
                end
                
                % Do final plot
                if(doplotfinal ==1)
                    subplot(2,2,3);
                    p.PlotArrayGeneric(obj.temp1, 'final plotting sensitivities after running. ')
                    drawnow
                end
                %             end
            end % end loading cases
            
            
            %                obj.temp1(ely,elx)  =    obj.temp1(ely,elx) /t2; % average the cases
            obj.temp1  =    obj.temp1 /t2; % average the cases
        end % end CalculateSensitiviesMesoStructureNoPeriodic
        
        
        
        % ------------------------------------------------------------
        % GetHomogenizedProperties
        % ------------------------------------------------------------
        %
        % For meso strudcture design HOmogenize the matrix and get K_h
        % Return the K_matrix that discribes this meso structure's
        % properties
        % ------------------------------------------------------------
        function macroElemProps = GetHomogenizedProperties(obj, config,homgSettings, matProp, loop,macroElemProps)
            
            [macroElemProps.D_homog]=FE_elasticV2_homgonization(obj, config, matProp);
            
        end % end CalculateSensitiviesMesoStructure
        
        
        
        
        %%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [dcn]=check(obj, nelx,nely,rmin,x,dc)
            dcn=zeros(nely,nelx);
            for i = 1:nelx
                for j = 1:nely
                    sum=0.0;
                    for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
                        for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
                            fac = rmin-sqrt((i-k)^2+(j-l)^2);
                            sum = sum+max(0,fac);
                            dcn(j,i) = dcn(j,i) + max(0,fac)*x(l,k)*dc(l,k);
                        end
                    end
                    dcn(j,i) = dcn(j,i)/(x(j,i)*sum);
                end
            end
        end
        
        %  -----------------------------------
        % Flip Exx and Eyy so that theta is always positive
        % %  -----------------------------------
        function [obj]=FlipOrientation(obj,config)
            
            for i = 1:config.nelx
                for j = 1:config.nely
                    if(obj.t(j,i)<0)
                        temp1=obj.Exx(j,i);
                        obj.Exx(j,i)=obj.Eyy(j,i);
                        obj.Eyy(j,i)=temp1;
                        obj.t(j,i)=obj.t(j,i)+pi/2;
                        
                        % Also flip the meso level values or else the
                        % optimization will be crazy!
%                         temp2 = obj.ExxSub(j,i);
%                         obj.ExxSub(j,i)=obj.EyySub(j,i);   % Sub system copies of design var
%                          obj.EyySub(j,i)=temp2;
%                         obj.thetaSub(j,i)= obj.thetaSub(j,i)+pi/2;
   
                        
                    end
                end
            end
        end
        
        
        
        % -----------------------------
        % GenerateStartingMesoDesign
        %
        % Geneate the x (density) values for the Meso design problem.
        % Several methods exist.
        % -----------------------------
        function [obj]= GenerateStartingMesoDesign(obj,mesoConfig,macroElementProperties)
            method =2;
            
            if(method ==1)
                obj.x(1:mesoConfig.nely,1:mesoConfig.nelx) = mesoConfig.totalVolume; % artificial density of the elements
                % method 1, randome values. Does not seem to be working well.
                obj.x(1:mesoConfig.nely,1:mesoConfig.nelx) = randi([0, round(mesoConfig.totalVolume*100)],mesoConfig.nely,mesoConfig.nelx)/100; % artificial density of the elements, can not be unifrom or else sensitivity will be 0 everywhere.
            elseif(method ==2)
                % method 2, box of empty in the middle.
                obj.x(1:mesoConfig.nely,1:mesoConfig.nelx) = ones(mesoConfig.nely,mesoConfig.nelx);
                midY = round(mesoConfig.nely/2);
                midX = round(mesoConfig.nelx/2);
                ratio = mesoConfig.nelx/mesoConfig.nely;
                vEmpty = mesoConfig.nelx*mesoConfig.nely-mesoConfig.totalVolume*mesoConfig.nelx*mesoConfig.nely;
                dimY = floor(sqrt(vEmpty/ratio));
                yStart = midY-floor(dimY/2);
                dimX =  floor(ratio*dimY);
                xStart = midX-floor(dimX/2);
                obj.x(yStart:yStart+dimY-1,xStart:xStart+dimX-1)= ones(dimY,dimX)*0.01;
                
                
            elseif(method ==3)
                % method 3, circle in the moddle
                obj.x(1:mesoConfig.nely,1:mesoConfig.nelx) = ones(mesoConfig.nely,mesoConfig.nelx);
                %            obj.x(1:mesoConfig.nely,1:mesoConfig.nelx) = randi([1, round(mesoConfig.totalVolume*100)],mesoConfig.nely,mesoConfig.nelx)/100; % artificial density of the elements, can not be unifrom or else sensitivity will be 0 everywhere.
                
                midY = round(mesoConfig.nely/2);
                midX = round(mesoConfig.nelx/2);
                mesoConfig.totalVolume=0.5;%*mesoConfig.nely*mesoConfig.nelx;
                radius = sqrt((mesoConfig.totalVolume*mesoConfig.nelx*mesoConfig.nely-mesoConfig.nelx*mesoConfig.nely)/(-pi));
                for i = 1:mesoConfig.nelx
                    for j = 1:mesoConfig.nely
                        %                 if sqrt((i-mesoConfig.nelx/2-0.5)^2+(j-mesoConfig.nely/2-0.5)*2) < min(mesoConfig.nelx,mesoConfig.nely)/3
                        %                     obj.x(j,i) = mesoConfig.totalVolume/2;
                        %                 end
                        d = sqrt((i-midX)^2+(j-midY)^2);
                        if(d<radius)
                            obj.x(j,i)= 0.01;
                        end
                    end
                end
                
            elseif(method ==4)
                % method 4, many circle holes
                obj.x(1:mesoConfig.nely,1:mesoConfig.nelx) = ones(mesoConfig.nely,mesoConfig.nelx);
                numHolesX = 7;
                numHolesY =7;
                totalHoles =numHolesX*numHolesY;
                XholeCenters = 1:mesoConfig.nelx/(numHolesX+1):mesoConfig.nelx;
                YholeCenters = 1:mesoConfig.nely/(numHolesY+1):mesoConfig.nely;
                
                radius = sqrt((mesoConfig.totalVolume*mesoConfig.nelx*mesoConfig.nely-mesoConfig.nelx*mesoConfig.nely)/(-pi*totalHoles));
                for i = 1:mesoConfig.nelx
                    for j = 1:mesoConfig.nely
                        %                 if sqrt((i-mesoConfig.nelx/2-0.5)^2+(j-mesoConfig.nely/2-0.5)*2) < min(mesoConfig.nelx,mesoConfig.nely)/3
                        %                     obj.x(j,i) = mesoConfig.totalVolume/2;
                        %                 end
                        for mm = 1:(1+numHolesX)
                            for nn = 1:(1+numHolesY)
                                midX= XholeCenters(mm);
                                midY= YholeCenters(nn);
                                
                                d = sqrt((i-midX)^2+(j-midY)^2);
                                if(d<radius)
                                    obj.x(j,i)=  0.01;
                                end
                            end
                        end
                    end
                end
            elseif(method ==5)
                % method 5, small hole in middle
                obj.x(1:mesoConfig.nely,1:mesoConfig.nelx) = ones(mesoConfig.nely,mesoConfig.nelx);
                midY = round(mesoConfig.nely/2);
                midX = round(mesoConfig.nelx/2);
                radius = 3;
                for i = 1:mesoConfig.nelx
                    for j = 1:mesoConfig.nely
                        %                 if sqrt((i-mesoConfig.nelx/2-0.5)^2+(j-mesoConfig.nely/2-0.5)*2) < min(mesoConfig.nelx,mesoConfig.nely)/3
                        %                     obj.x(j,i) = mesoConfig.totalVolume/2;
                        %                 end
                        d = sqrt((i-midX)^2+(j-midY)^2);
                        if(d<radius)
                            obj.x(j,i)=  0.01;
                        end
                        
                    end
                end
                
                % Randome circles
            elseif(method ==6)
                % method 4, many randome circle holes
                numHolesX =6;
                numHolesY =6;
                totalHoles =numHolesX*numHolesY;
                XholeCenters = randi([0, mesoConfig.nelx],1,numHolesX+1);
                YholeCenters =  randi([0, mesoConfig.nely],1,numHolesY+1);
                
                obj.x(1:mesoConfig.nely,1:mesoConfig.nelx) = ones(mesoConfig.nely,mesoConfig.nelx);
                
                radius = sqrt((mesoConfig.totalVolume*mesoConfig.nelx*mesoConfig.nely-mesoConfig.nelx*mesoConfig.nely)/(-pi*totalHoles));
                for i = 1:mesoConfig.nelx
                    for j = 1:mesoConfig.nely
                        %                 if sqrt((i-mesoConfig.nelx/2-0.5)^2+(j-mesoConfig.nely/2-0.5)*2) < min(mesoConfig.nelx,mesoConfig.nely)/3
                        %                     obj.x(j,i) = mesoConfig.totalVolume/2;
                        %                 end
                        for mm = 1:(1+numHolesX)
                            for nn = 1:(1+numHolesY)
                                midX= XholeCenters(mm);
                                midY= YholeCenters(nn);
                                
                                d = sqrt((i-midX)^2+(j-midY)^2);
                                if(d<radius)
                                    obj.x(j,i)=  0.01;
                                end
                            end
                            
                        end
                    end
                end
            end
        end
        
    end % End Methods
end % End Class DesignVars