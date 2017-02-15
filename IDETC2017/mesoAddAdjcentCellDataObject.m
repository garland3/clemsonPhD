classdef mesoAddAdjcentCellDataObject
    % stores information about the adjacent meso unit cells.
    %   Detailed explanation goes here
    
    properties        
        useAdjacent = 0; % 0 = no, 1 = true        
        % determines if to add the adjacent boundary of the cell in each of
        % the 4 directions.
        useUp=0;
        useDown=0;
        useRight=0;
        useLeft=0;
        
        % Row or column of the adacent cells
        upBoundary=0;
        downBoundary=0;
        rightBoundary=0;
        leftBoundary=0;
    end
    
    methods
        
        function DesignVars = AddAdjacentSensitivity(obj,settings, macroElemProps, DesignVars)
            
            minSensitivity = min(min(DesignVars.dc));
            
            if(obj.useUp ==1)
                index1 = (obj.upBoundary>=settings.voidMaterialDensityCutOff);
                index2 = (obj.upBoundary<settings.voidMaterialDensityCutOff);
                sensitivtyBoundary = zeros(size(obj.upBoundary));
                sensitivtyBoundary(index1)=minSensitivity;
               sensitivtyBoundary(index2)=0;
                DesignVars.dc(end,1:end) =  sensitivtyBoundary;
                DesignVars.dc(end-1,1:end) =  DesignVars.dc(end,1:end) ;
%             else
%                 DesignVars.dc(end,1:end) = 0;
%                 DesignVars.dc(end-1,1:end) =  0;
            end
            
            
            if(obj.useDown ==1)
                index1 = (obj.downBoundary>=settings.voidMaterialDensityCutOff);
                index2 = (obj.downBoundary<settings.voidMaterialDensityCutOff);
                  sensitivtyBoundary = zeros(size(obj.downBoundary));
                sensitivtyBoundary(index1)=minSensitivity;
               sensitivtyBoundary(index2)=0;               
                DesignVars.dc(1,1:end) =  sensitivtyBoundary;
                DesignVars.dc(2,1:end) =  DesignVars.dc(1,1:end);
%             else
%                 DesignVars.dc(1,1:end) =  0;
%                 DesignVars.dc(2,1:end) = 0;
            end
            
            if(obj.rightBoundary ==1)
                index1 = (obj.rightBoundary>=settings.voidMaterialDensityCutOff);
                index2 = (obj.rightBoundary<settings.voidMaterialDensityCutOff);
                sensitivtyBoundary = zeros(size(obj.rightBoundary));
                sensitivtyBoundary(index1)=minSensitivity;
               sensitivtyBoundary(index2)=0;               
                DesignVars.dc(1:end,end) =    sensitivtyBoundary;
                DesignVars.dc(1:end,end-1) =    DesignVars.dc(1:end,end);
%             else
%                 DesignVars.dc(1:end,end) = 0;
%                 DesignVars.dc(1:end,end-1) = 0;
            end
            
            if(obj.useLeft ==1)
                index1 = (obj.leftBoundary>=settings.voidMaterialDensityCutOff);
                index2 = (obj.leftBoundary<settings.voidMaterialDensityCutOff);
                sensitivtyBoundary = zeros(size(obj.leftBoundary));
                sensitivtyBoundary(index1)=minSensitivity;
               sensitivtyBoundary(index2)=0;   
                DesignVars.dc(1:end,1) =    sensitivtyBoundary;
                DesignVars.dc(1:end,2) =    sensitivtyBoundary;
%             else
%                 % the far left is fixed on the wall. 
%                 if(macroElemProps.xPosition ~=1)
%                     DesignVars.dc(1:end,1) = 0;
%                     DesignVars.dc(1:end,2) =  0;
%                 end
                
            end
            
        end
        
        
    end
    
end

