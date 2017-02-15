function  SingleMesoStuctureWrappe(settingscopy,ne,matProp)

folderNum = settingscopy.iterationNum;
macro_meso_iteration = settingscopy.macro_meso_iteration;
% GEt element->node mapping
outname = sprintf('./out%i/elementNodeMap%i.csv',folderNum,macro_meso_iteration);
IEN = csvread(outname);

%  read strains for all macro elements
% % Get displacement field
outname = sprintf('./out%i/displacement%i.csv',folderNum,macro_meso_iteration);
U =  csvread(outname);

% Get the density field, so we only need to use the strains with
% density above threshold
outname = sprintf('./out%i/densityfield%i.csv',folderNum,macro_meso_iteration);
x = csvread(outname);
% macroElementProps.density = x(macroElementProps.yPosition,macroElementProps.xPosition );

% Save element to XY position map (needed for x and w vars retrival)
outname = sprintf('./out%i/elementXYposition%i.csv',folderNum,macro_meso_iteration);
elementXYposition=csvread(outname);

[~, numMacroLoadingCases] = size(settingscopy.loadingCase);

macroElementProps = macroElementProp;
% find if the density is high enough for this element,
% if so, then add it as a loading case.
% count the loading cases.
loadcaseIndex = 1;
densityAverage=0;
elementUsedCount = 0;
for e = 1:ne
    results = elementXYposition(e,:);
    yPosition = results(1);
    xPosition = results(2);
    
    density = x(yPosition, xPosition);
    if(density<settingscopy.voidMaterialDensityCutOff)
        % if the density is less than 0.2, then don't even bother.
        continue;
    end
    densityAverage=densityAverage+density;
    elementUsedCount=elementUsedCount+1;
    
    nodes1=  IEN(e,:);
    macroElementProps.elementNodes=nodes1;
    xNodes = nodes1*2-1;
    yNodes = nodes1*2;
    dofNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
    
    % plan for multi-loading cases from the macro level.
    
    for loadcaseIndexMacro = 1:numMacroLoadingCases
        
        utemp = U(loadcaseIndexMacro,:);
        u_local =   utemp(dofNumbers);
        offsetX = u_local(7);
        offsetY = u_local(8);
        u_local([1 3 5 7]) = u_local([1 3 5 7])-offsetX;
        u_local([2 4 6 8]) = u_local([2 4 6 8])-offsetY;
        macroElementProps.disp(loadcaseIndex,:)  = u_local;
        loadcaseIndex=loadcaseIndex+1;
    end
    
end
% if(macro_meso_iteration>1)
outname = sprintf('./out%i/Dgiven_%i_forElement_%i.csv',folderNum,macro_meso_iteration,e);
D = csvread(outname);
macroElementProps.D_given =D;

densityAverage=densityAverage/elementUsedCount;


% mesoSettings.v1=matProp.CalculateDensityTargetforMeso(w,x,settingscopy);
% mesoSettings.v2=0;
% mesoSettings.totalVolume= mesoSettings.v1+0;

settings= settingscopy;
settings.loadingCase = zeros(1,loadcaseIndex-1);

settings.nelx = settings.nelxMeso;
settings.nely =settings.nelyMeso;


e = 1;
[DVMeso, ~] = GenerateDesignVarsForMesoProblem(settings,e);

macroElementProps.elementNumber = e;

% give periodic boundary condition.

DVMeso = DVMeso.CalcElementNodeMapmatrixWithPeriodicXandY(settings);
DVMeso =  DVMeso.CalcNodeLocationMeso(settings);

macroElementProps.material1Fraction=1;
[macroElementProps.K ,~,macroElementProps.B] = matProp.effectiveElasticKEmatrix( macroElementProps.material1Fraction, settings,[]);
macroElementProps.strain = macroElementProps.B* transpose(macroElementProps.disp); % transpose the disp to be vertical


 [D_homog,DVMeso,macroElementProps] = MesoStructureDesignV2(matProp,settings,DVMeso,macroElementProps,[]);

SaveMesoUnitCellDesignToCSV(DVMeso,macroElementProps,settings.iterationNum,settings.macro_meso_iteration,e,1);
% Homgenization(DVMeso, settings, matProp,macroElementProps,[]);
% disp(['Single Meso Design #: ' sprintf('%4i',e ) ' of ' sprintf('%4i',ne )]);
%  macroElementPropsParFor = GetMacroElementPropertiesFromCSV(settingscopy,e);
% scalePlot = 1;
% coord(:,1) = [0 1 1 0];
% coord(:,2)  = [0 0 1 1];
%
% % Check if void
% %         if(macroElementPropsParFor.density>settingscopy.voidMaterialDensityCutOff)
% if(macroElementPropsParFor.density>0.1)
%
%     % Only plot a few
%     settingscopy.mesoplotfrequency = 200;
%     if(mod(e,settingscopy.mesoplotfrequency) ==0)
%         settingscopy.doPlotAppliedStrain = 1;  plottingMesoDesign = 1;  plotting = 1; % this was for debugging % this was for debugging
%     else
%         settingscopy.doPlotAppliedStrain = 0; plottingMesoDesign = 0;    plotting = 0; % this was for debugging
%     end
%
%     % ------------------------------------------------------
%     % Single element per design var.
%     % ------------------------------------------------------
%     if(settingscopy.doUseMultiElePerDV==0)
%         if(plotting ==1)
%             figure(1)
%             [~, t2] = size(settingscopy.loadingCase);
%             for loadcaseIndex = 1:t2
%                 % utemp = U(loadcaseIndex,:);
%                 U2 = macroElementPropsParFor.disp(loadcaseIndex,:)*scalePlot;
%                 coordD = zeros(5,2);
%                 for temp = 1:4
%                     coordD(temp,1) =  coord(temp,1)+ U2(2*temp-1); % X value
%                     coordD(temp,2) =  coord(temp,2)+ U2(2*temp); % Y value
%                 end
%                 coord2 = coord;
%                 coordD(5,:) = coordD(1,:) ;
%                 coord2(5,:) = coord2(1,:);
%                 subplot(2,2*t2,2*loadcaseIndex-1);
%                 plot(coordD(:,1),coordD(:,2), '-b',coord2(:,1),coord2(:,2), '-g');
%                 axis([-0.3 1.3 -0.3 1.3])
%                 axis square
%             end
%         end
%
%     else
%         % ------------------------------------------------------
%         % Multiple element per design var.
%         % ------------------------------------------------------
%
%         if(e==-1)
%             macroElementPropsParFor.xDisplacements = [ 0 0 0 0 0 0  1 1 1]*0.1;
%             macroElementPropsParFor.yDisplacements = [ 0 0 0 0 0 0  1 1 1]*0.1;
%         end
%
%         if(plotting ==1)
%             figure(1)
%             [~, t2] = size(settingscopy.loadingCase);
%             for loadcaseIndex = 1:t2
%                 dx = macroElementPropsParFor.xDisplacements(loadcaseIndex,:)*scalePlot;
%                 dy = macroElementPropsParFor.yDisplacements(loadcaseIndex,:)*scalePlot;
%                 Xlocs = macroElementPropsParFor.mesoXnodelocations;
%                 Ylocs = macroElementPropsParFor.mesoYnodelocations;
%                 Xlocs = reshape(Xlocs',[],1);
%                 Ylocs = reshape(Ylocs',[],1);
%                 displacedX= Xlocs +dx';
%                 displacedY= Ylocs +dy';
%                 subplot(2,2*t2,2*loadcaseIndex-1);
%                 plot(Xlocs,Ylocs,'x');
%                 hold on
%                 plot(displacedX,displacedY,'o');
%                 hold off
%             end
%         end
%     end
%
%     [DVMeso, mesoSettings] = GenerateDesignVarsForMesoProblem(settingscopy,e);
%
%
%     % Set the target infill for the meso as the vol fraction of
%     %             mesoSettings.v1=0.5+(macroElementPropsParFor.material1Fraction*macroElementPropsParFor.density^settings.penal)/2;
%     w = macroElementPropsParFor.material1Fraction;
%     x = macroElementPropsParFor.density;
%     mesoSettings.v1=matProp.CalculateDensityTargetforMeso(w,x,settingscopy);
%     mesoSettings.v2=0;
%     mesoSettings.totalVolume= mesoSettings.v1+0;
%
%
%     mesoSettings.averageMultiElementStrain= settingscopy.averageMultiElementStrain;
%     mesoSettings.doPlotAppliedStrain=settingscopy.doPlotAppliedStrain;
% [D_homog,DVMeso,macroElementPropsParFor]= MesoStructureDesignV2(matProp,mesoSettings,DVMeso,macroElementPropsParFor,[]);
%     D_homog
%
%     newDesign = 1;
%
%
%     if(plottingMesoDesign ==1)
%         p = plotResults;
%         figure(1)
%         subplot(2,2,2);
%         outname = sprintf('meso structure for macro element %i density %f',e, mesoSettings.v1);
%         p.PlotArrayGeneric(DVMeso.x,outname);
%         %                 subplot(2,2,3);
%         %                 outname = sprintf('meso structure sensitivity %i density %f',e, mesoSettings.v1);
%         %                 p.PlotArrayGeneric(DVMeso.temp1,outname);
%         drawnow
%         nameGraph = sprintf('./out%i/elementpicture%i.png',settingscopy.iterationNum, e);
%         print(nameGraph,'-dpng')
%     end
% else
%     %D_homog_flat = zeros(1,9);
%     newDesign = 0; % false
%     DVMeso=[];
% end
%
% SaveMesoUnitCellDesignToCSV(DVMeso,macroElementPropsParFor,settingscopy.iterationNum,settingscopy.macro_meso_iteration,e,newDesign);
%
% % SavedDmatrix(e,:) = D_homog_flat;
%
% % write the density, volume fraction and topology fields to .csv files
% % make a list of elemenents that have material (ie, we don't need to
% % design the void regions)
% % loop over the elements with material.
% % 1. read the displacment field and vol frac, calculate the E_given
% % 2. Run the MesoStructureDesign
% % 3. Save the results of the meso-structure to a .csv file.