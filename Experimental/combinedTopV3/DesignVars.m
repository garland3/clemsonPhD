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
         
         % --------------------------
         % Support vars
         % -------------------------------
          xold; %                 
         temp1; % Sensitivity 1
         temp2; % Sensitivity 2
         dc; % Derivative of c (hence dc). C is the objective.    
         g1elastic; % Derivative of c with respect to a material change for the elastic
         g1heat; %  Derivative of c with respect to a material change for the heat
         IEN; % element to node map. Save this matrix, so it does not need to be recalculated every time.       
        XLocations; %=zeros(numNodesInRow,numNodesInColumn);
        YLocations; %=zeros(numNodesInRow,numNodesInColumn);
        globalPosition; %  = zeros(nn,2);
        
        
        NodeToXYArrayMap; % map of node numbers to their X,Y position in FEA arrays
         
     end
     
      methods
          % Constructur method
          function obj = DesignVars(settings)
              obj.CalcIENmatrix(settings);
              obj.CalcElementLocation(settings);
            
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
                    obj.IEN(count,:)=[rowMultiplier*elementsInRow+j,
                                  rowMultiplier*elementsInRow+j+1,
                                  (rowMultiplier +1)*elementsInRow+j+1,
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
          
          function [volume1, volume2] = CalculateVolumeFractions(obj)
              
              volume1 = 0;
              volume2 = 0;
              
              for i = 1:nelx
                for j = 1:nely
                    structureLocal = structure(j,i);
                    if(structureLocal == 0) % if void region
                        E_atElement(j,i) = E_empty;
                        K_atElement(i,j) = K_empty;
                        structGradArray(j,i) = Enylon-100;
                    else % if a filled region
                        volFraclocal = volFracArray(j,i);
                        volFracV1 = volFracV1 +volFraclocal; % sum up the total use of material 1 (PLA)
                        volFracV2 = volFracV2 + (1- volFraclocal); % sum up the total use of material 2 (Nylon)

                        K_atElement(i,j) = KheatPLA*volFraclocal+(1-volFraclocal)*KheatNylon;
                        E_atElement(j,i)= Epla*volFraclocal+(1-volFraclocal)*Enylon;  % simple mixture ratio
                        structGradArray(j,i) = E_atElement(j,i);
                    end
                end
            end

            % normalize the volume fraction by the number of elements
            volFracV1 = volFracV1/ne;
            volFracV2 = volFracV2/ne;
              
          end
          
      end
     
     
    
end