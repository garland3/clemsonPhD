function [DVmeso, ConfigMeso] = GenerateDesignVarsForMesoProblem(ConfigMeso,e,macroElementProperties)

% --------------------------------------
%  ConfigMeso
% --------------------------------------------

% target volumes of material 1 and 2
ConfigMeso.v1=macroElementProperties.targetDensity;
ConfigMeso.v2=0;
ConfigMeso.totalVolume= ConfigMeso.v1+0;

ConfigMeso.nelx = ConfigMeso.nelxMeso;
ConfigMeso.nely =ConfigMeso.nelyMeso;

% ---------------------------------
% Initialization of varriables
% ---------------------------------
DVmeso = DesignVars(ConfigMeso);

% Reuse the existing X matrix if it exists.
if(ConfigMeso.macro_meso_iteration>1)
    DesignNumber = ConfigMeso.macro_meso_iteration-1;
    folderNum = ConfigMeso.iterationNum;
    outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,DesignNumber,e);
    if exist(outname, 'file') == 2
        DVmeso.x = csvread(outname);
    else
       %designVars.x(1:ConfigMeso.nely,1:ConfigMeso.nelx) = ConfigMeso.totalVolume; % artificial density of the elements
         DVmeso.x(1:ConfigMeso.nely,1:ConfigMeso.nelx) = randi([0, ConfigMeso.totalVolume*100],ConfigMeso.nely,ConfigMeso.nelx)/100; % artificial density of the elements, can not be unifrom or else sensitivity will be 0 everywhere.
    end
else
DVmeso = DVmeso.PreCalculateXYmapToNodeNumber(ConfigMeso);
DVmeso.x(1:ConfigMeso.nely,1:ConfigMeso.nelx) =randi([0, 100],ConfigMeso.nely,ConfigMeso.nelx)/100*ConfigMeso.totalVolume; % artificial density of the elements, can not be unifrom or else sensitivity will be 0 everywhere.

    % designVars.x(1:ConfigMeso.nely,1:ConfigMeso.nelx) = ConfigMeso.totalVolume; % artificial density of the elements
    
    % method 1, randome values. Does not seem to be working well.
   
    
    % method 2, box of empty in the middle.
    %     designVars.x(1:ConfigMeso.nely,1:ConfigMeso.nelx) = ones(ConfigMeso.nely,ConfigMeso.nelx);
    %     midY = floor(ConfigMeso.nely/2); midX = floor(ConfigMeso.nelx/2);
    %     ratio = ConfigMeso.nelx/ConfigMeso.nely;
    %     vEmpty = ConfigMeso.nelx*ConfigMeso.nely-ConfigMeso.totalVolume*ConfigMeso.nelx*ConfigMeso.nely;
    %     dimY = floor(sqrt(vEmpty/ratio));
    %     yStart = midY-floor(dimY/2);
    %     dimX =  floor(ratio*dimY);
    %     xStart = midX-floor(dimX/2);
    %      designVars.x(yStart:yStart+dimY-1,xStart:xStart+dimX-1)= zeros(dimY,dimX);
    
    % method 3
    %     for i = 1:ConfigMeso.nelx
    %         for j = 1:ConfigMeso.nely
    %             if sqrt((i-ConfigMeso.nelx/2-0.5)^2+(j-ConfigMeso.nely/2-0.5)*2) < min(ConfigMeso.nelx,ConfigMeso.nely)/3
    %                 designVars.x(j,i) = ConfigMeso.totalVolume/2;
    %             end
    %         end
    %     end
    
end


% macro_meso_iteration=ConfigMeso.macro_meso_iteration;

% adjacent = mesoAddAdjcentCellDataObject;
% if(ConfigMeso.mesoAddAdjcentCellBoundaries==1)
%     
%    
%     % check if we should be using adjacent cells.
%      folderNum = ConfigMeso.iterationNum;
%   
%     if( ConfigMeso.useAjacentLocal==1 )
%         adjacent.useAdjacent = 1;
%         
%         outname = sprintf('./out%i/densityfield%i.csv',folderNum,macro_meso_iteration);
%         xMacro = csvread(outname);
%         
%         cutOff = ConfigMeso.noNewMesoDesignDensityCutOff;
%         DesignNumber = ConfigMeso.macro_meso_iteration;
%        
%         
%         % ------------------------------------
%         % check above
%         % ------------------------------------
%         if(macroElementProperties.yPosition <ConfigMeso.nely)
%             % Check the density of the adjacent up cell
%             yAdjacent=macroElementProperties.yPosition+1;
%             xAdjacent = macroElementProperties.xPosition;
%             if(xMacro(yAdjacent,xAdjacent)>cutOff)
%                 eAdjacent = macroElementProperties.elementNumber+ConfigMeso.nelx;
%                 % if dense enough then, get the lower row.
%                 outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,DesignNumber,eAdjacent);
%                 if exist(outname, 'file') == 2
%                     xDensityAdjacent = csvread(outname);
%                     adjacent.upBoundary=xDensityAdjacent(1,1:end);
%                     adjacent.useUp  =1;
%                 end
%             end
%         end
%         
%         % ------------------------------------
%         % check below
%         % ------------------------------------
%         if(macroElementProperties.yPosition >1)
%             % Check the density of the adjacent up cell
%             yAdjacent=macroElementProperties.yPosition-1;
%             xAdjacent = macroElementProperties.xPosition;
%             if(xMacro(yAdjacent,xAdjacent)>cutOff)
%                 eAdjacent = macroElementProperties.elementNumber-ConfigMeso.nelx;
%                 % if dense enough then, get the lower row.
%                 outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,DesignNumber,eAdjacent);
%                 if exist(outname, 'file') == 2
%                     xDensityAdjacent = csvread(outname);
%                     adjacent.downBoundary=xDensityAdjacent(end,1:end);
%                     adjacent.useDown =1;
%                 end
%             end
%             
%         end
%         
%         % ------------------------------------
%         % check right
%         % ------------------------------------
%         if(macroElementProperties.xPosition <ConfigMeso.nelx)
%             % Check the density of the adjacent up cell
%             yAdjacent=macroElementProperties.yPosition;
%             xAdjacent = macroElementProperties.xPosition+1;
%             if(xMacro(yAdjacent,xAdjacent)>cutOff)
%                 eAdjacent = macroElementProperties.elementNumber+1;
%                 % if dense enough then, get the lower row.
%                 outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,DesignNumber,eAdjacent);
%                 if exist(outname, 'file') == 2
%                     xDensityAdjacent = csvread(outname);
%                     adjacent.rightBoundary=xDensityAdjacent(1:end,1);
%                     adjacent.useRight =1;
%                 end
%             end
%         end
%         
%         % ------------------------------------
%         % check left
%         % ------------------------------------
%         if(macroElementProperties.xPosition >1)
%             % Check the density of the adjacent up cell
%             yAdjacent=macroElementProperties.yPosition;
%             xAdjacent = macroElementProperties.xPosition-1;
%             if(xMacro(yAdjacent,xAdjacent)>cutOff)
%                 eAdjacent = macroElementProperties.elementNumber-1;
%                 % if dense enough then, get the lower row.
%                 outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,DesignNumber,eAdjacent);
%                 if exist(outname, 'file') == 2
%                     xDensityAdjacent = csvread(outname);
%                     adjacent.leftBoundary=xDensityAdjacent(1:end,end);
%                     adjacent.useLeft =1;
%                 end
%             end
%         end
%     else
%         adjacent.useAdjacent = 0;
%     end  
% else
%      adjacent.useAdjacent = 0;
% end
%  designVars.mesoAddAdjcentCellDataObject = adjacent;
% if doing meso optimization, setup optimization configurations
%if ( ConfigMeso.mode == 4)

% designVars = designVars.CalcElementNodeMapmatrixWithPeriodicXandY(ConfigMeso);
% designVars =  designVars.CalcNodeLocationMeso(ConfigMeso);
%end





