classdef DesignVars
    
     properties
         x; % the "density" at each element
         xold; % 
         w; % the volume fraction at each element
         
         % Support vars to help in optimization
         temp1; % Sensitivity 1
         temp2; % Sensitivity 2
         dc; % Derivative of c (hence dc). C is the objective.
         
         IEN; % element to node map. Save this matrix, so it does not need to be recalculated every time. 
      
        XLocations; %=zeros(numNodesInRow,numNodesInColumn);
        YLocations; %=zeros(numNodesInRow,numNodesInColumn);
         
     end
     
      methods
          function obj = DesignVars(settings)
              obj.CalcIENmatrix(settings);
              obj.CalcElementLocation(settings);
             
              
              
          end
          
          function CalcElementLocation(obj,settings)
              
               numNodesInRow = settings.nelx + 1;
               numNodesInColumn = settings.nely + 1;
               obj.XLocations=zeros(numNodesInRow,numNodesInColumn);
               obj.YLocations=zeros(numNodesInRow,numNodesInColumn);
              
          end
          
          function CalcIENmatrix(obj,settings)
              
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
          
      end
     
     
    
end