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
         
         % --------------------------
         % Support vars
         % -------------------------------
          xold; %                 
         temp1; % Sensitivity 1
         temp2; % Sensitivity 2
         
         complianceSensitivity; % 
         dc; % Derivative of c (hence dc). C is the objective.    
         g1elastic; % Derivative of c with respect to a material change for the elastic
         g1heat; %  Derivative of c with respect to a material change for the heat
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
        
        
         
     end
     
      methods
          % Constructur method
          function obj = DesignVars(settings)
              obj.CalcIENmatrix(settings);
              obj.CalcElementLocation(settings);
              
              % Get the B matrix, the E,v,G do not matter and are not used
              % in the B calculation, so set them to 1. 
              E = 1; v= 1; G = 1;
              [~, ~, B_out] = elK_elastic(E,v, G);
              obj.B = B_out;
            
          end
          
          % Calcualte the Center of each element and put the information
          % into an array. Needed for the FEA
          % Calculate it here, so it only need to be calculated once. 
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
          
          function obj =  CalcIENmatrix(obj,settings)
               
            count = 1;
            elementsInRow = settings.nelx+1;
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
          
          function obj = CalculateSensitivies(obj, settings, matProp, loop)
              elementsInRow = settings.nelx+1;
              
              obj.xold = obj.x;
                % FE-ANALYSIS

                [obj.U_heatColumn]=temperatureFEA_V3(obj, settings, matProp,loop);   
                [obj.U, obj.maxF, obj.maxU]=FE_elasticV2(obj, settings, matProp);   
              

                % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
                    obj.c = 0.; % c is the objective. Total strain energy


                        for ely = 1:settings.nely
                            rowMultiplier = ely-1;
                            for elx = 1:settings.nelx

                                   nodes1=[rowMultiplier*elementsInRow+elx;
                                            rowMultiplier*elementsInRow+elx+1;
                                            (rowMultiplier +1)*elementsInRow+elx+1;
                                            (rowMultiplier +1)*elementsInRow+elx];

                                    % Get the element K matrix for this partiular element
                                    KE = matProp.effectiveElasticKEmatrix(  obj.w(ely,elx),settings);
                                    KEHeat = matProp.effectiveHeatKEmatrix(  obj.w(ely,elx), settings);
                                    Dmaterial = matProp.calculateEffectiveConstitutiveEquation( obj.w(ely,elx), settings);

                                    xNodes = nodes1*2-1;
                                    yNodes = nodes1*2;

                                 % NodeNumbers = union(xNodeNumbers,yNodeNumbers);
                                  NodeNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];

                                 % NodeNumbers = union(xNodeNumbers,yNodeNumbers);
                                  Ue = obj.U(NodeNumbers,:);
                                  U_heat = obj.U_heatColumn(nodes1,:);
                                  averageElementTemp = mean2(U_heat); % calculate the average temperature of the 4 nodes

                                  
                                   % obj.c =  obj.c + settings.w1*obj.x(ely,elx)^settings.penal*Ue'*KE*Ue;
                                   complianceTotal = settings.w1*obj.x(ely,elx)^settings.penal*Ue'*KE*Ue;
                                   
                                   % Find the elastic strain
                                   elasticStrain = obj.B*Ue;
                                   
                                   term1 = transpose(elasticStrain)*Dmaterial*elasticStrain*obj.x(ely,elx)^settings.penal;
                                   
                                   % Find the thermal strain 
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
                                   
                                   
                                   total =  settings.w1*(term1 + term2 + term3);
                                   [complianceTotal ,total]
                                   obj.c =  obj.c +total;
                                   obj.c =  obj.c + settings.w2*obj.x(ely,elx)^settings.penal*U_heat'*KEHeat*U_heat;     
                                   
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


                                % Topology sensitivies
                                %
                                % Calculate the sensitivity using the
                                % original method (with respec to min
                                % compliance)
                               obj.complianceSensitivity(ely,elx) = -settings.penal*obj.x(ely,elx)^(settings.penal-1)*Ue'*matProp.dKelastic*Ue;
                                
                               % calcualte the minimum strain energy
                               % sensitivity
                                obj.temp1(ely,elx) = -total;      
                               obj.temp1(ely,elx) = obj.complianceSensitivity(ely,elx);
                                
                                % calculate the minim temp sensitivity
                                 obj.temp2(ely,elx) = -settings.penal*obj.x(ely,elx)^(settings.penal-1)*U_heat'*KEHeat*U_heat;

                                 % Calculate the derivative with respect to a material
                                 % volume fraction composition change (not density change)
                                 obj.g1elastic(ely,elx) = obj.x(ely,elx)^(settings.penal)*Ue'*matProp.dKelastic*Ue;
                                 
                                 obj.g1heat(ely,elx) = obj.x(ely,elx)^(settings.penal)*U_heat'*matProp.dKheat*U_heat;
                            end
                        end
              
          end
          
      end
     
     
    
end