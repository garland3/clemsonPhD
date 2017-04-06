function [DVmeso,macroElemProps,finalDensity]=MesoStructureDesignV2(matProp,mesoConfig,DVmeso,macroElemProps,dcGiven)


doPlot =0; % For debugging allow plotting of some information.
if(mod(macroElemProps.elementNumber,mesoConfig.mesoplotfrequency) ==0)
    doPlot =1;
end
recvid = 0; % make a video?

% macroElemProps.disp= [ 0 0     0 0    0.1 0.1     0 0   ]

% Calcualte the strain, epsilon = B*d
% get the B matrix.
[~ ,~,macroElemProps.B] = matProp.effectiveElasticKEmatrix( macroElemProps.material1Fraction, mesoConfig,[]);
macroElemProps.strain = macroElemProps.B* transpose(macroElemProps.disp); % transpose the disp to be vertical

shearStrainSum= sum(macroElemProps.strain(3,:));

verboseOutput =0;
%  doPlot=1;


shearSign = 1;
if(shearStrainSum<0)
    shearSign=1;
else
    shearSign=-1;
end


if(doPlot ==1)
    p = plotResults;
    figure(1)
    subplot(1,1,1);
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

densityArray = 100;

old_muMatrix = ones(3,3);
penaltyValue=0.5;
macroElemProps.D_subSys = ones(3,3);

pstrain=ones(3,1)/3;
pstrain(3)=pstrain(3)*shearSign;
pstrainOld=ones(3,1);


counter = 1;
volumeUpdateInterval = 12;
mesoConfig.totalVolume=0.5; % start at 50% density and go from there. 

for mm = 1:mesoConfig.maxNumPseudoStrainLoop   
    
    temp=sum(abs(pstrain-pstrainOld));
    if(temp<mesoConfig.PseudoStrainEndCriteria && mm>1)
  
            
            temp2 = 'break becasue pseudo strain diff is small'       
             fprintf('Break Psudo Strain Loops value = [%f %f %f] PsuedoStrainLoop %i Diff %f\n',pstrain(1),pstrain(2),pstrain(3),mm,temp);
       
        break;
    end
    % ------------------------------------------------------
    % Generate the Design Vars for the Meso Optimization
    % ------------------------------------------------------
    DVmeso= DVmeso.GenerateStartingMesoDesign(mesoConfig,macroElemProps);
    
      pstrainOld=pstrain;
     oldX=  DVmeso.x+10;
    for mesoLoop = 1:mesoConfig.maxMesoLoops   
        macroElemProps.psuedoStrain=pstrain;
        [ DVmeso ,D_h, objective,muMatrix] = Homgenization(DVmeso, mesoConfig, matProp, macroElemProps,mesoLoop,old_muMatrix,penaltyValue);
        macroElemProps.D_subSys=D_h;
        
        
        % ------------------------------
        % Update the pseudo strain
        % ------------------------------
        pstrain=macroElemProps.psuedoStrain;
      
        diffSys=macroElemProps.D_sys./macroElemProps.D_subSys;
        damp =0.1;
        temp = [diffSys(1,1) ;diffSys(2,2); diffSys(3,3)];
        temp(3)=temp(3)*shearSign;
        temp = temp/(sum(abs(temp)));
        pstrain=pstrain+damp*temp;        
        
        pstrain=pstrain/(sum(sum(abs(pstrain))));        
        macroElemProps.psuedoStrain=pstrain;
        
        % ------------------------------
        % Update the volume constraint
        % ------------------------------
        counter=counter+1;
        if(mod(counter,volumeUpdateInterval)==1)
            termByTermDivision = macroElemProps.D_sys./macroElemProps.D_subSys;
            
            diff = -1+(termByTermDivision(1,1)+termByTermDivision(2,2))/2;
            damp=0.25;
            temp = diff*damp;
            maxChange = 0.1;
            temp = max(min(temp,maxChange),-maxChange);
            oldVolume =    mesoConfig.totalVolume;
            newVolume=mesoConfig.totalVolume+temp;
            minimumDensity=0.1;
            mesoConfig.totalVolume=min(max(newVolume,minimumDensity),1);
            ttttt=1;
        end       
                   
      
        % ------------------------------
        % Test for a converged design that is not changing. 
        % ------------------------------
           
        % test for convergeance of the x (density) values. Break out of
        % loop if they are not changing. 
        DiffX = oldX-DVmeso.x;
        sumOfDiffX = sum(sum(abs(DiffX)));
              if(verboseOutput==1)
          fprintf('%f \t%f \t %f and density %f, sumDiff %f\n',pstrain(1),pstrain(2),pstrain(3),      mesoConfig.totalVolume,sumOfDiffX);
              end
          changeLimit = mesoConfig.nelx*mesoConfig.nely*0.003;
        if(sumOfDiffX<changeLimit)
             if(verboseOutput==1)
              fprintf('Break in Meso design. X values are not changing. \n');
             end
           break; 
        end
        
         oldX=  DVmeso.x;        
    
        densityArray = [densityArray sum(sum(DVmeso.x))];       
     
        
        % FILTERING OF SENSITIVITIES
        [DVmeso.dc]   = DVmeso.check(mesoConfig.nelx,mesoConfig.nely,mesoConfig.rmin,DVmeso.x,DVmeso.dc);
        
        moveLimit=0.05;
        [DVmeso.x] = OC(mesoConfig.nelx,mesoConfig.nely,DVmeso.x,mesoConfig.totalVolume,DVmeso.dc, DVmeso, mesoConfig,moveLimit);
        
      
        
        
        if(doPlot ==1)
            figure(1)
            subplot(2,1,1)
            titleText = sprintf('meso design -> topology var Iter: %i',mesoLoop);
            p.PlotArrayGeneric(DVmeso.x,titleText); % plot the results.
          
            drawnow
        end
        
        
        if recvid==1
            drawnow
            F(vid) = getframe(figure(1)); % %Get frame of the topology in each iteration
            writeVideo(vidObj,F(vid)); %Save the topology in the video
            vid=vid+1;
        end
    end
end
% pstrain



if recvid==1
    close(vidObj);  %close video
end

% plot(densityArray)
% densityArray(end)

disp(['Meso Design #: ' sprintf('%4i',macroElemProps.elementNumber ) ' after '  sprintf('%4i',mesoLoop ) ' meso iterations: density= '  sprintf('%f',mesoConfig.totalVolume )]);

% --------------------------------------------
%    CALCULATE Effective constitutive matrix of the meso structure. This is
%    need for the macro optimization
% --------------------------------------------

% give periodic boundary condition.

% DVmeso = DVmeso.CalcElementNodeMapmatrixWithPeriodicXandY(mesoConfig);
% DVmeso =  DVmeso.CalcNodeLocationMeso(mesoConfig);

% macroElemProps = DVmeso.GetHomogenizedProperties(mesoConfig,mesoConfig, matProp, masterloop,macroElemProps);
macroElemProps.D_subSys=D_h;
finalDensity=  mesoConfig.totalVolume;
