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

addLineBetweeenDiverseMesoStructures = 0;

if(addLineBetweeenDiverseMesoStructures==1)
    folderName=sprintf('out%i',folderNum);
    outnameExx = sprintf('./%s/ExxSubSysValues%i.csv',folderName,macro_meso_iteration);
    outnameEyy = sprintf('./%s/EyySubSysValues%i.csv',folderName,macro_meso_iteration);
    outnametheta = sprintf('./%s/ThetaSubSysValues%i.csv',folderName,macro_meso_iteration);
    %  outnamerho = sprintf('./%s/densityUsedSubSysValues%i.csv',folderName,macro_meso_iteration);
    outnamerho = sprintf('./%s/densityUsedSubSysValues%i.csv',folderName,macro_meso_iteration);
    
    %------------------------
    % read the macro columns as well. This will help with future analysis
    % ----------------------------
    
    MesoExx=csvread(outnameExx);
    
    
    MesoEyy=csvread(outnameEyy);
    
    
    MesoTheta=csvread(outnametheta);
    
    
    MesoRho=csvread(outnamerho);
    
    [x,y] = meshgrid(1:config.nelx,1:config.nely);
    %     div = divergence(U,V)
    [px,py] = gradient(MesoExx);
    figure
    contour(x,y,MesoExx)
    hold on
    quiver(x,y,px,py)
    hold off
    nameGraph = sprintf('./GradientExxSub%i.png', config.macro_meso_iteration);
    print(nameGraph,'-dpng');
    
    [px,py] = gradient(MesoEyy);
    figure
    contour(x,y,MesoEyy)
    hold on
    quiver(x,y,px,py)
    hold off
    nameGraph = sprintf('./GradientEyySub%i.png', config.macro_meso_iteration);
    print(nameGraph,'-dpng');
    
    [px,py] = gradient(MesoTheta);
    figure
    contour(x,y,MesoTheta)
    hold on
    quiver(x,y,px,py)
    hold off
    nameGraph = sprintf('./GradientThetaSub%i.png', config.macro_meso_iteration);
    print(nameGraph,'-dpng');
    
end



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
%         fprintf('density of element %f\n',sum(sum(x)));
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
        fprintf('step 2 element %i of %i\n',e,ne);
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

if(config.multiscaleMethodCompare~=1)
    completeStruct( completeStruct>1)=1;
    
    completeStruct(completeStruct>config.voidMaterialDensityCutOff)=1;
    completeStruct(completeStruct<config.voidMaterialDensityCutOff)=0;
    % p.PlotArrayGenericWithBlueWhiteColors( completeStruct, 'Before')
    
    [completeStruct , numChanged] = CheckForConerElements(completeStruct, totalX,totalY, config.voidMaterialDensityCutOff);
    completeStruct(completeStruct<config.voidMaterialDensityCutOff)=0;
    
    completeStruct(completeStruct>config.voidMaterialDensityCutOff)=1;
    completeStruct(completeStruct<config.voidMaterialDensityCutOff)=0;
    p.PlotArrayGenericWithBlueWhiteColors( completeStruct, 'After Corner Check')
    
%     totalSize = size(completeStruct)
%     totalX
%     totalY
if 1==0
    sumBefore = sum(sum(completeStruct));
     [completeStruct ] = CheckRemoveOrphanedSegments(completeStruct, totalX,totalY, config.voidMaterialDensityCutOff);
    sumAfter = sum(sum(completeStruct));

    %   p.PlotArrayGenericWithBlueWhiteColors( completeStruct, 'After Isolated Element Check')
    
    fprintf('Change in density after removing elements is %f\n',sumBefore-sumAfter);
    end
end
% %end
density = sum(sum(completeStruct))/(totalX*totalY);
outname = sprintf('./mode90/Density%i.csv',config.macro_meso_iteration);
csvwrite(outname,density);
close all
figure(1)
plotname = sprintf('complete structure %i with density %f',config.macro_meso_iteration,density);
p.PlotArrayGenericWithBlueWhiteColors( completeStruct, plotname)

%     colorbar
freezeColors
nameGraph = sprintf('./completeStucture%f_macroIteration_%i.png', config.w1,config.macro_meso_iteration);
print(nameGraph,'-dpng', '-r1200')
if (config.generateCompleteStructureCSV==1)
    outname = sprintf('./completeStucture%f_macroIteration_%i.csv', config.w1,config.macro_meso_iteration);
    csvwrite(outname,completeStruct);
end