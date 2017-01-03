function [DVmeso,macroElemProps]=MesoStructureDesignV2(matProp,mesoConfig,DVmeso,macroElemProps,dcGiven)


doPlot =0; % For debugging allow plotting of some information.
 if(mod(macroElemProps.elementNumber,mesoConfig.mesoplotfrequency) ==0)
     doPlot =1;
 end
recvid = 0; % make a video?

% Calcualte the strain, epsilon = B*d
% get the B matrix.
[~ ,~,macroElemProps.B] = matProp.effectiveElasticKEmatrix( macroElemProps.material1Fraction, mesoConfig,[]);
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

DVmeso = DVmeso.CalcElementNodeMapmatrixWithPeriodicXandY(mesoConfig);
DVmeso =  DVmeso.CalcNodeLocationMeso(mesoConfig);


% --------------------------------------------
%    CALCULATE SENSITIVITY
%
%    UPDATE THE DESGIN OF THE UNIT CELL
% --------------------------------------------
% Loop Calculatint the sensitivity and changing the design var X
objectiveArray = 100;
densityArray = 100;
mesoConfig.penal = 3;
terminationMulti = 0.01;

% if(settings.singleMesoDesign~=1)
for mesoLoop = 1:100
    
    [ DVmeso ,D_h, objective] = Homgenization(DVmeso, mesoConfig, matProp, macroElemProps,mesoLoop);   
    macroElemProps.D_subSys=D_h;
    
%     objectiveArray = [objectiveArray objective];
%     densityArray = [densityArray sum(sum(DVmeso.x))];
    
    t = mesoConfig.terminationAverageCount;
    if(size(objectiveArray,2)>(t+2))    
        [avgObj] = FindAvergeChangeOfLastValue(objectiveArray, mesoConfig);
         [avgDensity] = FindAvergeChangeOfLastValue(densityArray, mesoConfig);
        if(abs(avgObj)<  mesoConfig.terminationCriteria*terminationMulti && avgDensity<  mesoConfig.terminationCriteria*terminationMulti )
            [avgObj avgDensity]
            break;
        end
    end
    
    
    
    
    
    % FILTERING OF SENSITIVITIES
    [DVmeso.dc]   = DVmeso.check(mesoConfig.nelx,mesoConfig.nely,mesoConfig.rmin,DVmeso.x,DVmeso.dc);
    
%     if(DVmeso.mesoAddAdjcentCellDataObject.useAdjacent ==1)
%        DVmeso= DVmeso.mesoAddAdjcentCellDataObject.AddAdjacentSensitivity(mesoConfig, macroElemProps, DVmeso);
%     end
    
    % DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
    [DVmeso.x] = OC(mesoConfig.nelx,mesoConfig.nely,DVmeso.x,mesoConfig.totalVolume,DVmeso.dc, DVmeso, mesoConfig);
    
    if(doPlot ==1)
        figure(1)
        subplot(2,2,3)
        titleText = sprintf('meso design -> topology var Iter: %i',mesoLoop);       
        p.PlotArrayGeneric(DVmeso.x,titleText); % plot the results.
        subplot(2,2,4)
        p.PlotArrayGeneric(DVmeso.dc,'meso design -> sensitivity'); % plot the results.
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
%     DVmeso.dcSum = zeros(settings.nelyMeso,settings.nelxMeso);
%     DVmeso.dc
%
%     for mesoLoop = 1:60
%         [ DVmeso ,D_h, objective] = Homgenization(DVmeso, mesoConfig, matProp, macroElemProps,mesoLoop);
%         change=objectiveold-objective;
%         objectiveold=objective;
%         relativeChange = change/objective;
%
%         if(abs(relativeChange)<  mesoConfig.terminationCriteria)
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

% DVmeso = DVmeso.CalcElementNodeMapmatrixWithPeriodicXandY(mesoConfig);
% DVmeso =  DVmeso.CalcNodeLocationMeso(mesoConfig);

% macroElemProps = DVmeso.GetHomogenizedProperties(mesoConfig,mesoConfig, matProp, masterloop,macroElemProps);
macroElemProps.D_subSys=D_h;
