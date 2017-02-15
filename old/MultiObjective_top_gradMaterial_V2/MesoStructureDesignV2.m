function [D_homog,designVarsMeso,macroElemProps]=MesoStructureDesignV2(matProp,mesoSettings,designVarsMeso,macroElemProps,dcGiven)



doPlot =0; % For debugging allow plotting of some information.
 if(mod(macroElemProps.elementNumber,mesoSettings.mesoplotfrequency) ==0)
     doPlot =1;
 end
recvid = 0; % make a video?

% Calcualte the strain, epsilon = B*d
% get the B matrix.
[macroElemProps.K ,~,macroElemProps.B] = matProp.effectiveElasticKEmatrix( macroElemProps.material1Fraction, mesoSettings,[]);
macroElemProps.strain = macroElemProps.B* transpose(macroElemProps.disp); % transpose the disp to be vertical

if(doPlot ==1)
    p = plotResults;
    figure(1)
    subplot(2,2,3);
end

if recvid==1
    videoOut = './resultsOuts.avi';
    vidObj = VideoWriter(videoOut);    %Prepare the new file for video
    vidObj.FrameRate = 5;
    vidObj.Quality = 100;
    open(vidObj);
    vid=1;
end

% give periodic boundary condition.

designVarsMeso = designVarsMeso.CalcElementNodeMapmatrixWithPeriodicXandY(mesoSettings);
designVarsMeso =  designVarsMeso.CalcNodeLocationMeso(mesoSettings);



% --------------------------------------------
%    CALCULATE SENSITIVITY
%
%    UPDATE THE DESGIN OF THE UNIT CELL
% --------------------------------------------
% Loop Calculatint the sensitivity and changing the design var X
objectiveArray = 100;
densityArray = 100;
mesoSettings.penal = 3;
terminationMulti = 0.01;

% if(settings.singleMesoDesign~=1)
for mesoLoop = 1:100
    
    [ designVarsMeso ,D_h, objective] = Homgenization(designVarsMeso, mesoSettings, matProp, macroElemProps,mesoLoop);   
    
    objectiveArray = [objectiveArray objective];
    densityArray = [densityArray sum(sum(designVarsMeso.x))];
    
    t = mesoSettings.terminationAverageCount;
    if(size(objectiveArray,2)>(t+2))    
        [avgObj] = FindAvergeChangeOfLastValue(objectiveArray, mesoSettings);
         [avgDensity] = FindAvergeChangeOfLastValue(densityArray, mesoSettings);
        if(abs(avgObj)<  mesoSettings.terminationCriteria*terminationMulti && avgDensity<  mesoSettings.terminationCriteria*terminationMulti )
            [avgObj avgDensity]
            break;
        end
    end
    
    
    
    
    
    % FILTERING OF SENSITIVITIES
    [designVarsMeso.dc]   = check(mesoSettings.nelx,mesoSettings.nely,mesoSettings.rmin,designVarsMeso.x,designVarsMeso.dc);
    
%     if(designVarsMeso.mesoAddAdjcentCellDataObject.useAdjacent ==1)
%        designVarsMeso= designVarsMeso.mesoAddAdjcentCellDataObject.AddAdjacentSensitivity(mesoSettings, macroElemProps, designVarsMeso);
%     end
    
    % DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
    [designVarsMeso.x] = OC(mesoSettings.nelx,mesoSettings.nely,designVarsMeso.x,mesoSettings.totalVolume,designVarsMeso.dc, designVarsMeso, mesoSettings);
    
    if(doPlot ==1)
        figure(1)
        subplot(2,2,3)
        p.PlotArrayGeneric(designVarsMeso.x,'meso design -> topology var'); % plot the results.
        subplot(2,2,4)
        p.PlotArrayGeneric(designVarsMeso.dc,'meso design -> sensitivity'); % plot the results.
        drawnow
    end
    
    
    if recvid==1
        drawnow
        F(vid) = getframe(figure(1)); % %Get frame of the topology in each iteration
        writeVideo(vidObj,F(vid)); %Save the topology in the video
        vid=vid+1;
    end
end
% else % settings.singleMesoDesign==1
%     % --------------------------------
%     % Single meso structure method
%     % --------------------------------
%     designVarsMeso.dcSum = zeros(settings.nelyMeso,settings.nelxMeso);
%     designVarsMeso.dc
%
%     for mesoLoop = 1:60
%         [ designVarsMeso ,D_h, objective] = Homgenization(designVarsMeso, mesoSettings, matProp, macroElemProps,mesoLoop);
%         change=objectiveold-objective;
%         objectiveold=objective;
%         relativeChange = change/objective;
%
%         if(abs(relativeChange)<  mesoSettings.terminationCriteria)
%             break;
%         end
%     end
%
% end

if recvid==1
    close(vidObj);  %close video
end

disp(['Meso Design #: ' sprintf('%4i',macroElemProps.elementNumber ) ' after '  sprintf('%4i',mesoLoop ) ' meso iterations']);

% --------------------------------------------
%    CALCULATE Effective constitutive matrix of the meso structure. This is
%    need for the macro optimization
% --------------------------------------------

% give periodic boundary condition.

% designVarsMeso = designVarsMeso.CalcElementNodeMapmatrixWithPeriodicXandY(mesoSettings);
% designVarsMeso =  designVarsMeso.CalcNodeLocationMeso(mesoSettings);

% macroElemProps = designVarsMeso.GetHomogenizedProperties(mesoSettings,mesoSettings, matProp, masterloop,macroElemProps);
macroElemProps.D_homog=D_h;
D_homog=D_h;