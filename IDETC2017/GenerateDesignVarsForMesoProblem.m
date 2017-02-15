function [DVmeso, mesoConfig] = GenerateDesignVarsForMesoProblem(mesoConfig,e,macroElementProperties)

% --------------------------------------
%  mesoConfig
% --------------------------------------------

% target volumes of material 1 and 2
mesoConfig.v1=macroElementProperties.targetDensity;
mesoConfig.v2=0;
mesoConfig.totalVolume= mesoConfig.v1+0;

mesoConfig.nelx = mesoConfig.nelxMeso;
mesoConfig.nely =mesoConfig.nelyMeso;

% ---------------------------------
% Initialization of varriables
% ---------------------------------
DVmeso = DesignVars(mesoConfig);

% % Reuse the existing X matrix if it exists.
% if(mesoConfig.macro_meso_iteration>1)
%     DesignNumber = mesoConfig.macro_meso_iteration-1;
%     folderNum = mesoConfig.iterationNum;
%     outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,DesignNumber,e);
%     if exist(outname, 'file') == 2
%         DVmeso.x = csvread(outname);
%     else
%         %DVmeso.x(1:mesoConfig.nely,1:mesoConfig.nelx) = mesoConfig.totalVolume; % artificial density of the elements
%         %  DVmeso.x(1:mesoConfig.nely,1:mesoConfig.nelx) = randi([0, mesoConfig.totalVolume*100],mesoConfig.nely,mesoConfig.nelx)/100; % artificial density of the elements, can not be unifrom or else sensitivity will be 0 everywhere.
%         %          DVmeso.x(1:mesoConfig.nely,1:mesoConfig.nelx) = zeros(mesoConfig.nely,mesoConfig.nelx);
%         
%         
%         
%     end
% else
    DVmeso = DVmeso.PreCalculateXYmapToNodeNumber(mesoConfig);
    % DVmeso.x(1:mesoConfig.nely,1:mesoConfig.nelx) =randi([0, 100],mesoConfig.nely,mesoConfig.nelx)/100*mesoConfig.totalVolume; % artificial density of the elements, can not be unifrom or else sensitivity will be 0 everywhere.
    
    method =3;
    if(method ==1)
        DVmeso.x(1:mesoConfig.nely,1:mesoConfig.nelx) = mesoConfig.totalVolume; % artificial density of the elements
        
        % method 1, randome values. Does not seem to be working well.
        DVmeso.x(1:mesoConfig.nely,1:mesoConfig.nelx) = randi([0, round(mesoConfig.totalVolume*100)],mesoConfig.nely,mesoConfig.nelx)/100; % artificial density of the elements, can not be unifrom or else sensitivity will be 0 everywhere.
        
    elseif(method ==2)
        % method 2, box of empty in the middle.
        DVmeso.x(1:mesoConfig.nely,1:mesoConfig.nelx) = ones(mesoConfig.nely,mesoConfig.nelx);
        midY = round(mesoConfig.nely/2);
        midX = round(mesoConfig.nelx/2);
        ratio = mesoConfig.nelx/mesoConfig.nely;
        vEmpty = mesoConfig.nelx*mesoConfig.nely-mesoConfig.totalVolume*mesoConfig.nelx*mesoConfig.nely;
        dimY = floor(sqrt(vEmpty/ratio));
        yStart = midY-floor(dimY/2);
        dimX =  floor(ratio*dimY);
        xStart = midX-floor(dimX/2);
        DVmeso.x(yStart:yStart+dimY-1,xStart:xStart+dimX-1)= ones(dimY,dimX);
        
        
    elseif(method ==3)
        % method 3, circle in the moddle
         DVmeso.x(1:mesoConfig.nely,1:mesoConfig.nelx) = ones(mesoConfig.nely,mesoConfig.nelx);
           DVmeso.x(1:mesoConfig.nely,1:mesoConfig.nelx) = randi([1, round(mesoConfig.totalVolume*100)],mesoConfig.nely,mesoConfig.nelx)/100; % artificial density of the elements, can not be unifrom or else sensitivity will be 0 everywhere.
      
          midY = round(mesoConfig.nely/2);
        midX = round(mesoConfig.nelx/2);
            radius = sqrt((mesoConfig.totalVolume*mesoConfig.nelx*mesoConfig.nely-mesoConfig.nelx*mesoConfig.nely)/(-pi));
                for i = 1:mesoConfig.nelx
                    for j = 1:mesoConfig.nely
        %                 if sqrt((i-mesoConfig.nelx/2-0.5)^2+(j-mesoConfig.nely/2-0.5)*2) < min(mesoConfig.nelx,mesoConfig.nely)/3
        %                     DVmeso.x(j,i) = mesoConfig.totalVolume/2;
        %                 end
                    d = sqrt((i-midX)^2+(j-midY)^2);
                    if(d<radius)
                          DVmeso.x(j,i)= 0.01;
                    end
        
        
                    end
                end
        
    elseif(method ==4)
        % method 4, many circle holes
        DVmeso.x(1:mesoConfig.nely,1:mesoConfig.nelx) = ones(mesoConfig.nely,mesoConfig.nelx);
        numHolesX = 7;
        numHolesY =7;
        totalHoles =numHolesX*numHolesY;
        XholeCenters = 1:mesoConfig.nelx/(numHolesX+1):mesoConfig.nelx;
        YholeCenters = 1:mesoConfig.nely/(numHolesY+1):mesoConfig.nely;
        
        
        
        radius = sqrt((mesoConfig.totalVolume*mesoConfig.nelx*mesoConfig.nely-mesoConfig.nelx*mesoConfig.nely)/(-pi*totalHoles));
        for i = 1:mesoConfig.nelx
            for j = 1:mesoConfig.nely
                %                 if sqrt((i-mesoConfig.nelx/2-0.5)^2+(j-mesoConfig.nely/2-0.5)*2) < min(mesoConfig.nelx,mesoConfig.nely)/3
                %                     DVmeso.x(j,i) = mesoConfig.totalVolume/2;
                %                 end
                for mm = 1:(1+numHolesX)
                    for nn = 1:(1+numHolesY)
                        midX= XholeCenters(mm);
                        midY= YholeCenters(nn);
                        
                        d = sqrt((i-midX)^2+(j-midY)^2);
                        if(d<radius)
                            DVmeso.x(j,i)=  0.01;
                        end
                    end
                    
                    
                end
            end
        end
      elseif(method ==5)
           % method 5, small hole in middle
            DVmeso.x(1:mesoConfig.nely,1:mesoConfig.nelx) = ones(mesoConfig.nely,mesoConfig.nelx);
          midY = round(mesoConfig.nely/2);
        midX = round(mesoConfig.nelx/2);
            radius = 3;
                for i = 1:mesoConfig.nelx
                    for j = 1:mesoConfig.nely
        %                 if sqrt((i-mesoConfig.nelx/2-0.5)^2+(j-mesoConfig.nely/2-0.5)*2) < min(mesoConfig.nelx,mesoConfig.nely)/3
        %                     DVmeso.x(j,i) = mesoConfig.totalVolume/2;
        %                 end
                    d = sqrt((i-midX)^2+(j-midY)^2);
                    if(d<radius)
                          DVmeso.x(j,i)=  0.01;
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
        
          DVmeso.x(1:mesoConfig.nely,1:mesoConfig.nelx) = ones(mesoConfig.nely,mesoConfig.nelx);
        
        radius = sqrt((mesoConfig.totalVolume*mesoConfig.nelx*mesoConfig.nely-mesoConfig.nelx*mesoConfig.nely)/(-pi*totalHoles));
        for i = 1:mesoConfig.nelx
            for j = 1:mesoConfig.nely
                %                 if sqrt((i-mesoConfig.nelx/2-0.5)^2+(j-mesoConfig.nely/2-0.5)*2) < min(mesoConfig.nelx,mesoConfig.nely)/3
                %                     DVmeso.x(j,i) = mesoConfig.totalVolume/2;
                %                 end
                for mm = 1:(1+numHolesX)
                    for nn = 1:(1+numHolesY)
                        midX= XholeCenters(mm);
                        midY= YholeCenters(nn);
                        
                        d = sqrt((i-midX)^2+(j-midY)^2);
                        if(d<radius)
                            DVmeso.x(j,i)=  0.01;
                        end
                    end
                    
                    
                end
            end
        end     
    end
% end

% macro_meso_iteration=mesoConfig.macro_meso_iteration;

% adjacent = mesoAddAdjcentCellDataObject;
% if(mesoConfig.mesoAddAdjcentCellBoundaries==1)
%
%
%     % check if we should be using adjacent cells.
%      folderNum = mesoConfig.iterationNum;
%
%     if( mesoConfig.useAjacentLocal==1 )
%         adjacent.useAdjacent = 1;
%
%         outname = sprintf('./out%i/densityfield%i.csv',folderNum,macro_meso_iteration);
%         xMacro = csvread(outname);
%
%         cutOff = mesoConfig.noNewMesoDesignDensityCutOff;
%         DesignNumber = mesoConfig.macro_meso_iteration;
%
%
%         % ------------------------------------
%         % check above
%         % ------------------------------------
%         if(macroElementProperties.yPosition <mesoConfig.nely)
%             % Check the density of the adjacent up cell
%             yAdjacent=macroElementProperties.yPosition+1;
%             xAdjacent = macroElementProperties.xPosition;
%             if(xMacro(yAdjacent,xAdjacent)>cutOff)
%                 eAdjacent = macroElementProperties.elementNumber+mesoConfig.nelx;
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
%                 eAdjacent = macroElementProperties.elementNumber-mesoConfig.nelx;
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
%         if(macroElementProperties.xPosition <mesoConfig.nelx)
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
%  DVmeso.mesoAddAdjcentCellDataObject = adjacent;
% if doing meso optimization, setup optimization configurations
%if ( mesoConfig.mode == 4)

% DVmeso = DVmeso.CalcElementNodeMapmatrixWithPeriodicXandY(mesoConfig);
% DVmeso =  DVmeso.CalcNodeLocationMeso(mesoConfig);
%end





