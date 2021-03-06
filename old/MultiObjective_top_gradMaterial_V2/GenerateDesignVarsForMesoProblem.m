function [designVarsMeso, meso_settings] = GenerateDesignVarsForMesoProblem(settings,e,macroElementProperties)

% --------------------------------------
% %% meso_settings
% --------------------------------------------
meso_settings = settings;
meso_settings.macro_meso_iteration = settings.macro_meso_iteration;
% copy a bunch of var set earlier
meso_settings.numXElmPerDV=settings.numXElmPerDV;
meso_settings.numYElmPerDV=    settings.numYElmPerDV;
meso_settings.doPlotAppliedStrain =   settings.doPlotAppliedStrain ;
meso_settings.loadingCase =   settings.loadingCase ;
meso_settings.doUseMultiElePerDV= settings.doUseMultiElePerDV; % Tell the meso settings about the mode.

meso_settings.elasticMaterialInterpMethod = 2; % Hashin�Shtrikam law (average of upper and lower boundary)
meso_settings.heatMaterialInterpMethod = 5; % Hashin�Shtrikam law (average of upper and lower boundary)

% target volumes of material 1 and 2
meso_settings.v1 = 0.2;
meso_settings.v2 = 0.2;


meso_settings.w1 = 1; % do not set to zero, instead set to 0.0001. Else we will get NA for temp2
meso_settings.iterationNum = 0;
meso_settings.doSaveDesignVarsToCSVFile = 0;
meso_settings.doPlotFinal = 0;
meso_settings.terminationCriteria =settings.terminationCriteria; % 10%

% if meso structure designing, then make a smaller initial mesh

meso_settings.nelx = settings.nelxMeso;
meso_settings.nely =settings.nelyMeso;

meso_settings= meso_settings.UpdateVolTargetsAndObjectiveWeights();
%meso_settings


% ---------------------------------
% Initialization of varriables
% ---------------------------------
designVars = DesignVars(meso_settings);

% Reuse the existing X matrix if it exists.
if(settings.macro_meso_iteration>1)
    DesignNumber = settings.macro_meso_iteration-1;
    folderNum = settings.iterationNum;
    outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,DesignNumber,e);
    if exist(outname, 'file') == 2
        designVars.x = csvread(outname);
    else
        designVars.x(1:meso_settings.nely,1:meso_settings.nelx) = meso_settings.totalVolume; % artificial density of the elements
    end
    
else
    % designVars.x(1:meso_settings.nely,1:meso_settings.nelx) = meso_settings.totalVolume; % artificial density of the elements
    
    % method 1, randome values. Does not seem to be working well.
    designVars.x(1:meso_settings.nely,1:meso_settings.nelx) = randi([0, meso_settings.totalVolume*100],meso_settings.nely,meso_settings.nelx)/100; % artificial density of the elements, can not be unifrom or else sensitivity will be 0 everywhere.
    
    % method 2, box of empty in the middle.
    %     designVars.x(1:meso_settings.nely,1:meso_settings.nelx) = ones(meso_settings.nely,meso_settings.nelx);
    %     midY = floor(meso_settings.nely/2); midX = floor(meso_settings.nelx/2);
    %     ratio = meso_settings.nelx/meso_settings.nely;
    %     vEmpty = meso_settings.nelx*meso_settings.nely-meso_settings.totalVolume*meso_settings.nelx*meso_settings.nely;
    %     dimY = floor(sqrt(vEmpty/ratio));
    %     yStart = midY-floor(dimY/2);
    %     dimX =  floor(ratio*dimY);
    %     xStart = midX-floor(dimX/2);
    %      designVars.x(yStart:yStart+dimY-1,xStart:xStart+dimX-1)= zeros(dimY,dimX);
    
    % method 3
    %     for i = 1:meso_settings.nelx
    %         for j = 1:meso_settings.nely
    %             if sqrt((i-meso_settings.nelx/2-0.5)^2+(j-meso_settings.nely/2-0.5)*2) < min(meso_settings.nelx,meso_settings.nely)/3
    %                 designVars.x(j,i) = meso_settings.totalVolume/2;
    %             end
    %         end
    %     end
    
end
designVars.w(1:meso_settings.nely,1:meso_settings.nelx)  = 1; % actual volume fraction composition of each element
% fractionCurrent_V1Local =1;

% designVars.temp1(1:meso_settings.nely,1:meso_settings.nelx) = 0;
% designVars.temp2(1:meso_settings.nely,1:meso_settings.nelx) = 0;
% designVars.complianceSensitivity(1:meso_settings.nely,1:meso_settings.nelx) = 0;
% designVars.totalStress(1:meso_settings.nely,1:meso_settings.nelx) = 0;

% designVars.g1elastic(1:meso_settings.nely,1:meso_settings.nelx) = 0;
% designVars.g1heat(1:meso_settings.nely,1:meso_settings.nelx) = 0;

% designVars = designVars.CalcIENmatrix(meso_settings);
% designVars =  designVars.CalcNodeLocation(meso_settings);
designVars = designVars.PreCalculateXYmapToNodeNumber(meso_settings);
macro_meso_iteration=settings.macro_meso_iteration;

 adjacent = mesoAddAdjcentCellDataObject;
if(meso_settings.mesoAddAdjcentCellBoundaries==1)
    
   
    % check if we should be using adjacent cells.
     folderNum = settings.iterationNum;
  
    if( settings.useAjacentLocal==1 )
        adjacent.useAdjacent = 1;
        
        outname = sprintf('./out%i/densityfield%i.csv',folderNum,macro_meso_iteration);
        xMacro = csvread(outname);
        
        cutOff = settings.noNewMesoDesignDensityCutOff;
        DesignNumber = settings.macro_meso_iteration;
       
        
        % ------------------------------------
        % check above
        % ------------------------------------
        if(macroElementProperties.yPosition <settings.nely)
            % Check the density of the adjacent up cell
            yAdjacent=macroElementProperties.yPosition+1;
            xAdjacent = macroElementProperties.xPosition;
            if(xMacro(yAdjacent,xAdjacent)>cutOff)
                eAdjacent = macroElementProperties.elementNumber+settings.nelx;
                % if dense enough then, get the lower row.
                outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,DesignNumber,eAdjacent);
                if exist(outname, 'file') == 2
                    xDensityAdjacent = csvread(outname);
                    adjacent.upBoundary=xDensityAdjacent(1,1:end);
                    adjacent.useUp  =1;
                end
            end
        end
        
        % ------------------------------------
        % check below
        % ------------------------------------
        if(macroElementProperties.yPosition >1)
            % Check the density of the adjacent up cell
            yAdjacent=macroElementProperties.yPosition-1;
            xAdjacent = macroElementProperties.xPosition;
            if(xMacro(yAdjacent,xAdjacent)>cutOff)
                eAdjacent = macroElementProperties.elementNumber-settings.nelx;
                % if dense enough then, get the lower row.
                outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,DesignNumber,eAdjacent);
                if exist(outname, 'file') == 2
                    xDensityAdjacent = csvread(outname);
                    adjacent.downBoundary=xDensityAdjacent(end,1:end);
                    adjacent.useDown =1;
                end
            end
            
        end
        
        % ------------------------------------
        % check right
        % ------------------------------------
        if(macroElementProperties.xPosition <settings.nelx)
            % Check the density of the adjacent up cell
            yAdjacent=macroElementProperties.yPosition;
            xAdjacent = macroElementProperties.xPosition+1;
            if(xMacro(yAdjacent,xAdjacent)>cutOff)
                eAdjacent = macroElementProperties.elementNumber+1;
                % if dense enough then, get the lower row.
                outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,DesignNumber,eAdjacent);
                if exist(outname, 'file') == 2
                    xDensityAdjacent = csvread(outname);
                    adjacent.rightBoundary=xDensityAdjacent(1:end,1);
                    adjacent.useRight =1;
                end
            end
        end
        
        % ------------------------------------
        % check left
        % ------------------------------------
        if(macroElementProperties.xPosition >1)
            % Check the density of the adjacent up cell
            yAdjacent=macroElementProperties.yPosition;
            xAdjacent = macroElementProperties.xPosition-1;
            if(xMacro(yAdjacent,xAdjacent)>cutOff)
                eAdjacent = macroElementProperties.elementNumber-1;
                % if dense enough then, get the lower row.
                outname = sprintf('./out%i/densityfield%iforElement%i.csv',folderNum,DesignNumber,eAdjacent);
                if exist(outname, 'file') == 2
                    xDensityAdjacent = csvread(outname);
                    adjacent.leftBoundary=xDensityAdjacent(1:end,end);
                    adjacent.useLeft =1;
                end
            end
        end
    else
        adjacent.useAdjacent = 0;
    end  
else
     adjacent.useAdjacent = 0;
end
 designVars.mesoAddAdjcentCellDataObject = adjacent;
% if doing meso optimization, setup optimization configurations
%if ( meso_settings.mode == 4)

% designVars = designVars.CalcElementNodeMapmatrixWithPeriodicXandY(meso_settings);
% designVars =  designVars.CalcNodeLocationMeso(meso_settings);
%end



designVarsMeso=designVars;

