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
terminationMulti = 0.01;

% if(settings.singleMesoDesign~=1)
oldTerm1Max=1;
oldTerm2Max=1000;
lambdaMatrix = ones(3,3);
macroElemProps.D_subSys = ones(3,3);

macroElemProps.D_sys

%   mesoConfig.totalVolume=0.95;
%   mesoConfig.rmin=5;
% mesoConfig.maxMesoLoops=30;
% for matVolLoop =1:20
    for mesoLoop = 1:mesoConfig.maxMesoLoops
        
        [ DVmeso ,D_h, objective,term1Max,term2Max,lambdaMatrix,dx3] = Homgenization(DVmeso, mesoConfig, matProp, macroElemProps,mesoLoop,oldTerm1Max,oldTerm2Max,lambdaMatrix);
        macroElemProps.D_subSys=D_h;
        oldTerm1Max=term1Max;
        oldTerm2Max=term2Max;
        
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
      
     
         [DVmeso.x] = OC_meso(mesoConfig.nelx,mesoConfig.nely,DVmeso.x,mesoConfig.totalVolume,DVmeso.dc, DVmeso, mesoConfig);
         
%         move = 0.1;
%           [dx3]   = DVmeso.check(mesoConfig.nelx,mesoConfig.nely,mesoConfig.rmin,DVmeso.x,dx3);
%         x =DVmeso.x; 
%         x_iplus1 = x-DVmeso.dc./dx3;
%         xnew = max(0.01,max(x-move,min(1.,min(x+move,x_iplus1))));   
%         DVmeso.x = xnew;


%           NE = mesoConfig.nelx*mesoConfig.nely;
%          lambda = 1./DVmeso.dc;
%         move = 0.1;
%           [dx3]   = DVmeso.check(mesoConfig.nelx,mesoConfig.nely,mesoConfig.rmin,DVmeso.x,dx3);
%         x =DVmeso.x; 
%         x_iplus1 = x- ((NE + lambda.*DVmeso.dc)./(lambda.*dx3));
%         xnew = max(0.01,max(x-move,min(1.,min(x+move,x_iplus1))));   
%         DVmeso.x = xnew;
     
        
        if(doPlot ==1)
            figure(1)
            subplot(2,2,3)
            titleText = sprintf('meso design -> topology var Iter: %i',mesoLoop);
            p.PlotArrayGeneric(DVmeso.x,titleText); % plot the results.
            subplot(2,2,4)
            p.PlotArrayGeneric(DVmeso.dc,'meso design -> sensitivity'); % plot the results.
            
            subplot(2,2,1)
            p.PlotArrayGeneric(dx3,'meso design -> dx3'); % plot the results.
            drawnow
        end
        
        
        if recvid==1
            drawnow
            F(vid) = getframe(figure(1)); % %Get frame of the topology in each iteration
            writeVideo(vidObj,F(vid)); %Save the topology in the video
            vid=vid+1;
        end
    end
    
%     Diff_Sys_Sub =  (macroElemProps.D_sys-macroElemProps.D_subSys);
%     diff = sum(sum(Diff_Sys_Sub));
%     
%     terms=macroElemProps.D_sys./macroElemProps.D_subSys;
%     avg4terms = (terms(1,1)+terms(1,2)+terms(2,2)+terms(3,3))/4;
%     
%     if(diff>100)
%         % needs more material
%         macroElemProps.targetDensity=min(macroElemProps.targetDensity+0.1*diff/matProp.E_material1,1);
% %          macroElemProps.targetDensity=min(macroElemProps.targetDensity*avg4terms,1);
%     elseif(diff<100)
%         
%         macroElemProps.targetDensity=max(macroElemProps.targetDensity-0.1*abs(diff)/matProp.E_material1,0.05);
% %              macroElemProps.targetDensity=max(macroElemProps.targetDensity*avg4terms,0.05);
%     else
%         break;
%     end
%     mesoConfig.maxMesoLoops=mesoConfig.maxMesoLoops+2;
%     mesoConfig.totalVolume=macroElemProps.targetDensity;
%     density =   mesoConfig.totalVolume
%   
%    % method 3, circle in the moddle
%          DVmeso.x(1:mesoConfig.nely,1:mesoConfig.nelx) = ones(mesoConfig.nely,mesoConfig.nelx);
%            DVmeso.x(1:mesoConfig.nely,1:mesoConfig.nelx) = randi([0, round(mesoConfig.totalVolume*100)],mesoConfig.nely,mesoConfig.nelx)/100; % artificial density of the elements, can not be unifrom or else sensitivity will be 0 everywhere.
%       
%           midY = round(mesoConfig.nely/2);
%         midX = round(mesoConfig.nelx/2);
%             radius = sqrt((mesoConfig.totalVolume*mesoConfig.nelx*mesoConfig.nely-mesoConfig.nelx*mesoConfig.nely)/(-pi));
%                 for i = 1:mesoConfig.nelx
%                     for j = 1:mesoConfig.nely
%         %                 if sqrt((i-mesoConfig.nelx/2-0.5)^2+(j-mesoConfig.nely/2-0.5)*2) < min(mesoConfig.nelx,mesoConfig.nely)/3
%         %                     DVmeso.x(j,i) = mesoConfig.totalVolume/2;
%         %                 end
%                     d = sqrt((i-midX)^2+(j-midY)^2);
%                     if(d<radius)
%                           DVmeso.x(j,i)= 0;
%                     end
%         
%         
%                     end
%                 end
% end
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
