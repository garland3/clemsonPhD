function  MesoDesignWrapper(config,e,ne,matProp)


macroElementProperties = GetMacroElementPropertiesFromCSV(config,e);
disp(['Meso Design #: ' sprintf('%4i',e ) ' of ' sprintf('%4i',ne ) ...
    ' position X = '  sprintf('%4i',macroElementProperties.xPos) ' Y = ' sprintf('%4i',macroElementProperties.yPos) ...
    ' Target Density =  ' sprintf('%4i', macroElementProperties.targetDensity) ...
    ' MesoMacro Iteration =  ' sprintf('%4i', config.macro_meso_iteration) ]);
scalePlot = 1;
coord(:,1) = [0 1 1 0];
coord(:,2)  = [0 0 1 1];

% ----------------------------------------------
% Check Density
% Check if void we actually need to make a new design.
% -----------------------------------------------
if(macroElementProperties.densitySIMP>config.noNewMesoDesignDensityCutOff)
   
    % ------------------------------------------------------
    % Plot (only a few)
    % ------------------------------------------------------
     if(mod(e,config.mesoplotfrequency) ==0)
        config.doPlotAppliedStrain = 1;  plottingMesoDesign = 1;  plotting = 1; % this was for debugging % this was for debugging
    else
        config.doPlotAppliedStrain = 0; plottingMesoDesign = 0;    plotting = 0; % this was for debugging
    end
    
    if(plotting ==1)
        figure(1)
        [~, t2] = size(config.loadingCase);
        for loadcaseIndex = 1:t2
            % utemp = U(loadcaseIndex,:);
            U2 = macroElementProperties.disp(loadcaseIndex,:)*scalePlot;
            coordD = zeros(5,2);
            for temp = 1:4
                coordD(temp,1) =  coord(temp,1)+ U2(2*temp-1); % X value
                coordD(temp,2) =  coord(temp,2)+ U2(2*temp); % Y value
            end
            coord2 = coord;
            coordD(5,:) = coordD(1,:) ;
            coord2(5,:) = coord2(1,:);
            subplot(2,2*t2,2*loadcaseIndex-1);
            plot(coordD(:,1),coordD(:,2), '-b',coord2(:,1),coord2(:,2), '-g');
            axis([-0.3 1.3 -0.3 1.3])
            axis square
        end
    end
    
    % ------------------------------------------------------
    % Generate the Design Vars for the Meso Optimization
    % ------------------------------------------------------
    [DVmeso, configMeso] = GenerateDesignVarsForMesoProblem(config,e,macroElementProperties);
    
    % ----------------------------------
    % Actual Meso Design. 
    % ----------------------------------
    [DVmeso,macroElementProperties]= MesoStructureDesignV2(matProp,configMeso,DVmeso,macroElementProperties,[]);
    
   
    macroElementProperties.D_subSys  % Show the new D found by the design. 
    Diff_Sys_Sub =  (macroElementProperties.D_subSys- macroElementProperties.D_sys);
    Diff_Sys_Sub
    
    newDesign = 1;
    if(plottingMesoDesign ==1)
        p = plotResults;
        figure(1)
        subplot(2,2,2);
        outname = sprintf('meso structure for macro element %i density %f',e, configMeso.v1);
        p.PlotArrayGeneric(DVmeso.x,outname);
        %                 subplot(2,2,3);
        %                 outname = sprintf('meso structure sensitivity %i density %f',e, configMeso.v1);
        %                 p.PlotArrayGeneric(DVmeso.temp1,outname);
        drawnow
        nameGraph = sprintf('./out%i/elementpicture%i.png',config.iterationNum, e);
        print(nameGraph,'-dpng')
    end
else
    %D_homog_flat = zeros(1,9);
    newDesign = 0; % false
    DVmeso=[];
end

SaveMesoUnitCellDesignToCSV(DVmeso,macroElementProperties,config.iterationNum,config.macro_meso_iteration,e,newDesign);
