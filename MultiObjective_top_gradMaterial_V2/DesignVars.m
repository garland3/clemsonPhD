classdef DesignVars
    % Design varriables and support temp varriables are stored in this class.
    
    properties
        % --------------------------------
        % Design Var Arrays
        % --------------------------------
        x; % the "density" at each element
        w; % the volume fraction at each element
        
        % Optimization vars
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
        temp1; % Sensitivity 1
        temp2; % Sensitivity 2
        
        %complianceSensitivity; %
        totalStress;
        dc; % Derivative of c (hence dc). C is the objective.
        g1elastic; % Derivative of c with respect to a material change for the elastic
        g1heat; %  Derivative of cHeat with respect to a material change for the heat
        IEN; % element to node map. Save this matrix, so it does not need to be recalculated every time.
        XLocations; %=zeros(numNodesInRow,numNodesInColumn);
        YLocations; %=zeros(numNodesInRow,numNodesInColumn);
        globalPosition; %  = zeros(nn,2);
        NodeToXYArrayMap; % map of node numbers to their X,Y position in FEA arrays
        
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
    end
    
    methods
        % Constructur method
        function obj = DesignVars(settings)     
            
            % Get the B matrix, the E,v,G do not matter and are not used
            % in the B calculation, so set them to 1.
            E = 1; v= 1; G = 1; strain = [];
            [~, ~, B_out] = elK_elastic(E,v, G,strain,[]);
            obj.B = B_out;
            
        end
        
        
        
        % Calcualte the Center of each element and put the information
        % into an array. Needed for the FEA
        % Calculate it here, so it only need to be calculated once.
        %
        % Think these are actually node locations
        function obj = CalcElementLocation(obj,settings)
            nn = (settings.nelx+1)*(settings.nely+1); % number of nodes
            
            numNodesInRow = settings.nelx + 1;
            numNodesInColumn = settings.nely + 1;
            obj.XLocations=zeros(numNodesInRow,numNodesInColumn);
            obj.YLocations=zeros(numNodesInRow,numNodesInColumn);
            
            obj.globalPosition = zeros(nn,2);
            count = 1;
            for i = 1:numNodesInColumn  % y
                for j= 1:numNodesInRow % x
                    obj.globalPosition(count,:) = [j-1 i-1];
                    count = count+1;
                    obj.XLocations(j,i) = j-1;
                    obj.YLocations(j,i) = i-1;
                end
            end
        end
        
        
        % Calcualte the Center of each element and put the information
        % into an array. Needed for the FEA
        % Calculate it here, so it only need to be calculated once.
        function obj = CalcNodeLocationMeso(obj,settings)
            nn = (settings.nelx)*(settings.nely); % number of nodes
            
            numNodesInRow = settings.nelx ;
            numNodesInColumn = settings.nely ;
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
        
        function obj = CalcNodeLocationMeso_Tile(obj,settings)
            nn = (obj.nelxTile)*(obj.nelyTile); % number of nodes same as number of elements since it wraps
            %             nn = (settings.nelx)*(settings.nely); % number of nodes
            
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
        function obj =  CalcIENmatrix(obj,settings)
            
            count = 1;
            elementsInRow = settings.nelx+1; % think this is actually "nodes in a row. "
            nn = (settings.nelx+1)*(settings.nely+1); % number of nodes
            obj.IEN = zeros(nn,4);
            % Each row, so nely # of row
            for i = 1:settings.nely
                rowMultiplier = i-1;
                % Each column, so nelx # of row
                for j= 1:settings.nelx
                    obj.IEN(count,:)=[rowMultiplier*elementsInRow+j, ...
                        rowMultiplier*elementsInRow+j+1, ...
                        (rowMultiplier +1)*elementsInRow+j+1,...
                        (rowMultiplier +1)*elementsInRow+j];
                    count = count+1;
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
        function obj =  CalcElementNodeMapmatrixWithPeriodicXandY(obj,settings)
            
            count = 1;
            elementsInRow = settings.nelx;
            nn = (settings.nelx)*(settings.nely); % number of nodes same as number of elements since it wraps
            obj.IEN = zeros(nn,4);
            % Each row, so nely # of row
            for i = 1:settings.nely
                rowMultiplier = i-1;
                % Each column, so nelx # of row
                for j= 1:settings.nelx
                    
                    % normal case
                    if(j ~= settings.nelx && i ~= settings.nely )
                        obj.IEN(count,:)=[rowMultiplier*elementsInRow+j, ...
                            rowMultiplier*elementsInRow+j+1, ...
                            (rowMultiplier +1)*elementsInRow+j+1,...
                            (rowMultiplier +1)*elementsInRow+j];
                    elseif(j == settings.nelx && i ~= settings.nely )
                        % On the right side of the mesh case
                        
                        % loop back arround case case
                        obj.IEN(count,:)=[rowMultiplier*elementsInRow+j, ...
                            rowMultiplier*elementsInRow+1, ... % note the difference
                            (rowMultiplier +1)*elementsInRow+1,... % note the difference
                            (rowMultiplier +1)*elementsInRow+j];
                    elseif(j ~= settings.nelx && i == settings.nely )
                        % On the top side of the mesh case
                        
                        obj.IEN(count,:)=[rowMultiplier*elementsInRow+j, ...
                            rowMultiplier*elementsInRow+j+1, ...
                            0+j+1,... % note the difference
                            0+j]; % note the difference
                        
                    elseif(j == settings.nelx && i == settings.nely )
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
        function obj =  CalcElementNodeMapmatrixWithPeriodicXandY_Tile(obj,settings)
            
            count = 1;
            elementsInRow=obj.nelxTile;% elementsInRow = settings.nelx;
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
        function number = GetNodeNumberGivenXY(settings, x,y)
            numNodesInRow = settings.nelx+1;
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
        
        % Pre Calculate a map of the XY array coordinates for each node
        % number.
        function obj= PreCalculateXYmapToNodeNumber(obj ,settings)
            nn = (settings.nelx+1)*(settings.nely+1); % number of nodes
            obj.NodeToXYArrayMap = zeros(nn,2);
            count = 1;
            for i = 1:settings.nely
                for j= 1:settings.nelx
                    obj.NodeToXYArrayMap(count,:) = [i,j];
                    count = count+1;
                end
            end
        end
        
        % Calculates the volume
        function [volume1, volume2] = CalculateVolumeFractions(obj, settings)
            
            volume1 = 0;
            volume2 = 0;
            
            for i = 1:settings.nelx
                for j = 1:settings.nely
                    x_local = obj.x(j,i);
                    if(x_local <= settings.voidMaterialDensityCutOff) % if void region
                        % E_atElement(j,i) = E_empty;
                        % K_atElement(i,j) = K_empty;
                        % structGradArray(j,i) = Enylon-100;
                    else % if a filled region
                        volFraclocal = obj.w(j,i);
                        volume1 = volume1 +volFraclocal; % sum up the total use of material 1
                        volume2 = volume2 + (1- volFraclocal); % sum up the total use of material 2
                        
                        % K_atElement(i,j) = KheatPLA*volFraclocal+(1-volFraclocal)*KheatNylon;
                        % E_atElement(j,i)= Epla*volFraclocal+(1-volFraclocal)*Enylon;  % simple mixture ratio
                        % structGradArray(j,i) = E_atElement(j,i);
                    end
                end
            end
            
            % normalize the volume fraction by the number of elements
            ne = settings.nelx*settings.nely;
            volume1 = volume1/ne;
            volume2 = volume2/ne;
            
        end
        %%
        
        function obj = CalculateSensitivies(obj, settings, matProp, loop,macro_meso_iteration)
            Dgiven = [];
            elementsInRow = settings.nelx+1;
            %             obj.xold = obj.x;
            
            % FE-ANALYSIS
            [obj.U_heatColumn]=temperatureFEA_V3(obj, settings, matProp,loop);
            [obj.U, obj.maxF, obj.maxU]=FE_elasticV2(obj, settings, matProp);
            
            % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
            obj.c = 0.; % c is the objective. Total strain energy
            obj.cCompliance = 0;
            obj.cHeat = 0;
            
            count =1;
            
            for ely = 1:settings.nely
                rowMultiplier = ely-1;
                for elx = 1:settings.nelx
                    nodes1=[rowMultiplier*elementsInRow+elx;
                        rowMultiplier*elementsInRow+elx+1;
                        (rowMultiplier +1)*elementsInRow+elx+1;
                        (rowMultiplier +1)*elementsInRow+elx];
                    
                    xNodes = nodes1*2-1;
                    yNodes = nodes1*2;
                    NodeNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
                    
                    Ue = obj.U(NodeNumbers,:);
                    U_heat = obj.U_heatColumn(nodes1,:);
                    averageElementTemp = mean2(U_heat); % calculate the average temperature of the 4 nodes
                    
                    % Get the element K matrix for this partiular element
                    if(macro_meso_iteration>1)
                        e = count;
                        Dgiven =matProp.GetSavedDMatrix(e);
                    end
                    KE = matProp.effectiveElasticKEmatrix(  obj.w(ely,elx),settings,Dgiven);
                    KEHeat = matProp.effectiveHeatKEmatrix(  obj.w(ely,elx), settings);
                    Dmaterial = matProp.calculateEffectiveConstitutiveEquation( obj.w(ely,elx), settings,Dgiven);
                    
                    % Find the elastic strain
                    elasticStrain = obj.B*Ue;
                    term1 = transpose(elasticStrain)*Dmaterial*elasticStrain*obj.x(ely,elx)^settings.penal;
                    
                    % Find the thermal strain and add this in to the
                    if(settings.addThermalExpansion ==1)
                        alpha = matProp.effectiveThermalExpansionCoefficient(  obj.w(ely,elx));
                        deltaTemp = averageElementTemp- settings.referenceTemperature;
                        thermalStrain = alpha*deltaTemp*[1; 1; 0];
                        term2 = transpose(thermalStrain)*Dmaterial*thermalStrain*obj.x(ely,elx)^settings.penal;
                        term3 = -2*transpose(thermalStrain)*Dmaterial*elasticStrain*obj.x(ely,elx)^settings.penal;
                    else
                        term2=0;
                        term3=0;
                        thermalStrain=0;
                    end
                    
                    % Sum the elastic compliance terms.   % Compliance,
                    total = (term1 + term2 + term3);                   
                    
                    obj.cCompliance = obj.cCompliance + obj.x(ely,elx)^settings.penal*Ue'*KE*Ue;
                    obj.cHeat =   obj.cHeat           + obj.x(ely,elx)^settings.penal*U_heat'*KEHeat*U_heat;
                    
                    % Derivative of  D
                    % (constitutive matrix) with respect to
                    % a density change
                    dD =  settings.penal*obj.x(ely,elx)^(settings.penal-1)*Dmaterial;
                    dTerm1 = transpose(elasticStrain)*dD*elasticStrain;
                    
                    if(settings.addThermalExpansion ==1)
                        dTerm2 =transpose(thermalStrain)*dD*thermalStrain;
                        dTerm3 = -2*transpose(thermalStrain)*dD*elasticStrain;
                    else
                        dTerm2=0;
                        dTerm3=0;
                    end
                    
                    if (settings.doPlotStress == 1)
                        totalStrain=thermalStrain+elasticStrain;
                        totalSressLocal=Dmaterial*totalStrain*obj.x(ely,elx)^settings.penal;
                        vonM = sqrt(totalSressLocal(1)^2  +   totalSressLocal(2)^2 -  totalSressLocal(1)*totalSressLocal(2)  +   3*(totalSressLocal(3))^2);
                        obj.totalStress(ely,elx)=vonM;
                    end
                    
                    % Topology sensitivies
                    %
                    % Calculate the sensitivity using the
                    % original method (with respec to min
                    % compliance)
                    % obj.complianceSensitivity(ely,elx) = -settings.penal*obj.x(ely,elx)^(settings.penal-1)*Ue'*matProp.dKelastic*Ue;
                    
                    % calcualte the minimum strain energy
                    % sensitivity
                    
                    obj.temp1(ely,elx) = -total;
                    %obj.temp1(ely,elx) = obj.complianceSensitivity(ely,elx);
                    
                    % calculate the minim temp sensitivity
                    obj.temp2(ely,elx) = -settings.penal*obj.x(ely,elx)^(settings.penal-1)*U_heat'*KEHeat*U_heat;                    
                    
                    % Calculate the derivative with respect to a material
                    % volume fraction composition change (not density change)
                    totalMaterialD = dTerm1+dTerm2+dTerm3;
                    obj.g1elastic(ely,elx) =totalMaterialD;
                    %  obj.g1elastic(ely,elx) = obj.x(ely,elx)^(settings.penal)*Ue'*matProp.dKelastic*Ue;
                    
                    obj.g1heat(ely,elx) = obj.x(ely,elx)^(settings.penal)*U_heat'*matProp.dKheat*U_heat;
                    count=count+1;
                end % end, loop over nelx
            end % end, for loop over nely
        end % End Function, CalculateSenstivities
        
        
        % ------------------------------------------------------------
        %
        % CalculateSensitiviesMesoStructure
        %
        % ------------------------------------------------------------
        function obj = CalculateSensitiviesMesoStructure(obj, settings, matProp, loop,macroElemProps, U)
            
            doplot = 0;
            if(doplot ==1)
                p = plotResults;
            end
            
            % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
            obj.c = 0.; % c is the objective. Total strain energy
            obj.cCompliance = 0;
            obj.cHeat = 0;
            
            ne = settings.nelx*settings.nely; % number of elements
            for e = 1:ne
                
                % loop over local node numbers to get their node global node numbers
                nodes1 = obj.IEN(e,:);
                [elx,ely]= obj.GivenNodeNumberGetXY(e);
                
                
                xNodes = nodes1*2-1;
                yNodes = nodes1*2;
                dofNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
                
                Ue = U(dofNumbers);
                % U_heat = obj.U_heatColumn(nodes1,:);
                %averageElementTemp = mean2(U_heat); % calculate the average temperature of the 4 nodes
                
                % Get the element K matrix for this partiular element
                KE = matProp.effectiveElasticKEmatrix(  obj.w(ely,elx),settings,[]);
             
                % KEHeat = matProp.effectiveHeatKEmatrix(  obj.w(ely,elx), settings);
                % Dmaterial = matProp.calculateEffectiveConstitutiveEquation( obj.w(ely,elx), settings);
%                 settings.nelx
                % Find the elastic strain
                % elasticStrain = obj.B*Ue;
                term1 = transpose(Ue)*KE*Ue*obj.x(ely,elx)^settings.penal;
                %  term2 = 0;
                % term3= 0;
                
                % Sum the elastic compliance terms.
                % total = (term1 + term2 + term3);
                %
                
                obj.temp1(ely,elx) = -term1;
                
                if(doplot ==1)
                    if(mod(e,10) ==0)
                        p.PlotArrayGeneric(obj.temp1, 'plotting sensitivities while running. ')
                        drawnow
                    end
                end
                % calculate the minim temp sensitivity
                % obj.temp2(ely,elx) = -settings.penal*obj.x(ely,elx)^(settings.penal-1)*U_heat'*KEHeat*U_heat;
            end
            
            % Do final plot
            if(doplot ==1)
                p.PlotArrayGeneric(obj.temp1, 'final plotting sensitivities while running. ')
                    drawnow
             end
            %             end
        end % end CalculateSensitiviesMesoStructure
        
        % --------------------------------------------
        % The matrix given is 9 repeating unit cells. Only use the middle
        % one
        % same thing as above, but use a tiled version.
        % --------------------------------------------
        function obj = CalculateSensitiviesMesoStructure_Tile(obj, settings, matProp, loop,macroElemProps, U)
            doplot = settings.plotSensitivityWhilerunning;
%             doplot = 1;
            if(doplot ==1)
                p = plotResults;
            end
            % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
            obj.c = 0.; % c is the objective. Total strain energy
            obj.cCompliance = 0;
            obj.cHeat = 0;
            
%             ne = obj.nelyTile *obj.nelxTile ; % number of elements
%             elementsPerTile = settings.nely*settings.nelx;
%             start = (settings.sensitivityTile-1)*elementsPerTile;
%             endElement = (settings.sensitivityTile)*elementsPerTile;
            %             for e = start:endElement
            
            
            ne = settings.nelx*settings.nely; % number of elements
            
               eleInRow = settings.nelx*settings.numTilesX;
            
            % loop over the single unit cell ,and get sensitivity 
            for e = 1:ne
                
%                 columnsUp = settings.nely+floor(e/settings.nely);
                columnsUp = settings.nely+floor((e-1)/settings.nely);
                offset = eleInRow*columnsUp  +  settings.nelx;
                
                rowsOver = mod(e-1,settings.nelx)+1;
                tiledElementNum = rowsOver+offset;
                eTile=tiledElementNum;
                
                % loop over local node numbers to get their node global node numbers
                nodes1 = obj.IENTile(eTile,:);
                results = obj.globalPositionTile(eTile,:);
                elx=results(1); ely= results(2);
                
                %                 [elx,ely]= obj.GivenNodeNumberGetXY(e);
                
                xNodes = nodes1*2-1;
                yNodes = nodes1*2;
                dofNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
                Ue = U(dofNumbers); % Ue = obj.U(dofNumbers,:);
                % U_heat = obj.U_heatColumn(nodes1,:);
                %averageElementTemp = mean2(U_heat); % calculate the average temperature of the 4 nodes
                
                % Get the element K matrix for this partiular element
                KE = matProp.effectiveElasticKEmatrix(  obj.wTile(ely,elx),settings,[]);
                % KEHeat = matProp.effectiveHeatKEmatrix(  obj.w(ely,elx), settings);
                % Dmaterial = matProp.calculateEffectiveConstitutiveEquation( obj.w(ely,elx), settings);
                
                % Find the elastic strain
                % elasticStrain = obj.B*Ue;
                term1 = transpose(Ue)*KE*Ue*obj.xTile(ely,elx)^settings.penal;
                %  term2 = 0;
                % term3= 0;
                % Sum the elastic compliance terms.
                % total = (term1 + term2 + term3);
                
              result = obj.globalPosition(e,:);
              elxSingle=result(1)+1;    elySingle=result(2)+1;
%                 [elxSingle,elySingle]= obj.GivenNodeNumberGetXY(e);
                obj.temp1(elySingle,elxSingle) = -term1;
                
                if(doplot ==1)
                    p.PlotArrayGeneric(obj.temp1, 'plotting sensitivities while running. ')
                    drawnow
                end
                % calculate the minim temp sensitivity
                % obj.temp2(ely,elx) = -settings.penal*obj.x(ely,elx)^(settings.penal-1)*U_heat'*KEHeat*U_heat;
            end
%                         p = plotResults;
%                         p.PlotArrayGeneric( obj.temp1,'sensitivity');
        end % end CalculateSensitiviesMesoStructure
        
        % ------------------------------------------------------------
        % GetHomogenizedProperties
        % ------------------------------------------------------------
        %
        % For meso strudcture design HOmogenize the matrix and get K_h
        % Return the K_matrix that discribes this meso structure's
        % properties
        % ------------------------------------------------------------
        function macroElemProps = GetHomogenizedProperties(obj, settings,homgSettings, matProp, loop,macroElemProps)
            % Test 3 loading cases. XX, YY, XY
            %              loadstrain = [1 0 0;
            %                            0 1 0;
            %                            0 0 1]
            
            % load  = 1:3
            % for ll = load
            % I call the constiut
            [macroElemProps.D_homog]=FE_elasticV2_homgonization(obj, settings, matProp);
            %end
        end % end CalculateSensitiviesMesoStructure
        
    end % End Methods
end % End Class DesignVars