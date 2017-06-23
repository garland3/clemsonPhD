% --------------------------------------------------------
% --------------------------------------------------------
%
%
%           Genereate the complete macro-meso design!
%
%
% --------------------------------------------------------
% --------------------------------------------------------
function []= GenerateCompleteStructureV2Improved(config)


postProcess = 1;

close all
p = plotResults;


temp = config.mesoAddAdjcentCellBoundaries;
config.mesoAddAdjcentCellBoundaries=0;
macroElementProps = macroElementProp;
macroElementProps.targetDensity=0.5; % make up some value for now. 
[DVMeso, mesoconfig] = GenerateDesignVarsForMesoProblem(config,1,macroElementProps);
config.mesoAddAdjcentCellBoundaries=temp;

mesoconfig.doUseMultiElePerDV =config.doUseMultiElePerDV;
numTilesX=config.numTilesX;
numTilesY = config.numTilesY;

% Generate huge area
totalX=config.nelx*mesoconfig.nelx*numTilesX
totalY=config.nely*mesoconfig.nely*numTilesY

completeStruct = zeros(totalY,totalX);
ne = config.nelx*config.nely; % number of elements


%--------------------------------------------
% Get the density field
%--------------------------------------------
macro_meso_iteration = config.macro_meso_iteration;
%macroElementProps = macroElementProp;
% macroElementProps.elementNumber = e;
folderNum = config.iterationNum;
% GET the saved element to XY position map (needed for x and w vars retrival)
outname = sprintf('./out%i/elementXYposition%i.csv',folderNum,macro_meso_iteration);
elementXYposition=csvread(outname);
% Get the density field
outname = sprintf('./out%i/SIMPdensityfield%i.csv',folderNum,macro_meso_iteration);
xxx = csvread(outname);




for e = 1:ne
    fprintf('element %i of %i\n',e,ne);
    macroElementProps.elementNumber=e;
    results = elementXYposition(macroElementProps.elementNumber,:);
    macroElementProps.yPos = results(1);
    macroElementProps.xPos = results(2);
    macroElementProps.densitySIMP = xxx(macroElementProps.yPos,macroElementProps.xPos );
    
    % Check if void
    if(macroElementProps.densitySIMP>config.voidMaterialDensityCutOff || config.multiscaleMethodCompare)
        x=GetMesoUnitCellDesignFromCSV(config,e);
        DVMeso.x = x;
        % -------------------------------------
        % No, post processing
        % -------------------------------------
        if(postProcess~=1)
            yShift = (macroElementProps.yPos-1)*mesoconfig.nely*numTilesY+1;
            xShift = (macroElementProps.xPos-1)*mesoconfig.nelx*numTilesX+1;
            DVMeso=TileMesoStructure(mesoconfig, DVMeso);
            completeStruct(yShift:(yShift+mesoconfig.nely*numTilesY-1),xShift:(xShift+mesoconfig.nelx*numTilesX-1))=DVMeso.xTile;
        else
            % -------------------------------------
            % Yes, With, post processing
            % -------------------------------------
            step = 1;
            completeStruct= TileMesoStructureV2(mesoconfig,config, DVMeso,macroElementProps,xxx,completeStruct,step);
            
        end
        
        
    end
end

if(postProcess==1)
    for e = 1:ne
        fprintf('step 2element %i of %i\n',e,ne);
        macroElementProps.elementNumber=e;
        results = elementXYposition(macroElementProps.elementNumber,:);
        macroElementProps.yPos = results(1);
        macroElementProps.xPos = results(2);
        macroElementProps.densitySIMP = xxx(macroElementProps.yPos,macroElementProps.xPos );
        
        % Check if void
        %if(macroElementProps.densitySIMP>config.voidMaterialDensityCutOff)
        step = 2;
        completeStruct= TileMesoStructureV2(mesoconfig,config, DVMeso,macroElementProps,xxx,completeStruct,step);
        %end
    end
end



%if(config.multiscaleMethodCompare~=1)
    % set the max value to be 1
%      completeStruct=    imgaussfilt(completeStruct,2 );
    completeStruct( completeStruct>1)=1;
   
    completeStruct(completeStruct>config.voidMaterialDensityCutOff)=1;
    completeStruct(completeStruct<config.voidMaterialDensityCutOff)=0;
    
     [completeStruct , numChanged] = CheckForConerElements(completeStruct, totalX,totalY, config.voidMaterialDensityCutOff);
      completeStruct(completeStruct<config.voidMaterialDensityCutOff)=0;
      
      [completeStruct ] = CheckRemoveOrphanedSegments(completeStruct, totalX,totalY, config.voidMaterialDensityCutOff);
%end
density = sum(sum(completeStruct))/(totalX*totalY);
 outname = sprintf('./mode90/Density%i.csv',config.macro_meso_iteration);
    csvwrite(outname,density);

figure(1)
plotname = sprintf('complete structure %i with density %f',config.macro_meso_iteration,density);
p.PlotArrayGeneric( completeStruct, plotname)
rgbSteps = 100;  caxis([0,1]);
map = colormap; % current colormap
middPoint = floor(rgbSteps/4);
map(1:middPoint,:) = [ones(middPoint,1),ones(middPoint,1),ones(middPoint,1)];
for zz =    middPoint:rgbSteps
    map(zz,:) = [0,               1- zz/rgbSteps, 0.5];
end
colormap(map)
%     colorbar
freezeColors
nameGraph = sprintf('./completeStucture%f_macroIteration_%i.png', config.w1,config.macro_meso_iteration);
print(nameGraph,'-dpng', '-r1200')
if (config.generateCompleteStructureCSV==1)
    outname = sprintf('./completeStucture%f_macroIteration_%i.csv', config.w1,config.macro_meso_iteration);
    csvwrite(outname,completeStruct);
end