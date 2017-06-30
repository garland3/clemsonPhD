function  MesoDesignWrapper(config,e,ne,matProp)


fprintf('PreRead Macro State for Meso Element %i\n',e)
macroElementProperties = GetMacroElementPropertiesFromCSV(config,e);
if(macroElementProperties.targetDensity<0)
    return;
end
disp(['Meso Design #: ' sprintf('%4i',macroElementProperties.elementNumber  ) ' of ' sprintf('%4i',ne ) ...
    ' position X = '  sprintf('%4i',macroElementProperties.xPos) ' Y = ' sprintf('%4i',macroElementProperties.yPos) ...
    ' MesoMacro Iteration =  ' sprintf('%4i', config.macro_meso_iteration) ]);
% scalePlot = 1;
% coord(:,1) = [0 1 1 0];
% coord(:,2)  = [0 0 1 1];
 if(config.strainAndTargetTest==1 || config.UseLookUpTableForPsuedoStrain==1)
     p=macroElementProperties.psuedoStrain;
     fprintf('Meso design , pseudo strain %f %f %f density target %f\n',p(1),p(2),p(3),macroElementProperties.targetDensity);
%           fprintf('Meso design as ANN training data, pseudo strain %f %f %f density target %f\n',p(1),p(2),p(3),macroElementProperties.targetDensity);
 end

% ----------------------------------------------
% Check Density
% Check if void we actually need to make a new design.
% -----------------------------------------------
if(macroElementProperties.densitySIMP>config.noNewMesoDesignDensityCutOff || config.multiscaleMethodCompare==1)
    
    % ------------------------------------------------------
    % Plot (only a few)
    % ------------------------------------------------------
    if(mod(e,config.mesoplotfrequency) ==0)
        config.doPlotAppliedStrain = 1;  plottingMesoDesign = 1;  plotting = 1; % this was for debugging % this was for debugging
    else
        config.doPlotAppliedStrain = 0; plottingMesoDesign = 0;    plotting = 0; % this was for debugging
    end
    
    %     if(plotting ==1)
    %         figure(1)
    %         [~, t2] = size(config.loadingCase);
    %         for loadcaseIndex = 1:t2
    %             % utemp = U(loadcaseIndex,:);
    %             U2 = macroElementProperties.disp(loadcaseIndex,:)*scalePlot;
    %             coordD = zeros(5,2);
    %             for temp = 1:4
    %                 coordD(temp,1) =  coord(temp,1)+ U2(2*temp-1); % X value
    %                 coordD(temp,2) =  coord(temp,2)+ U2(2*temp); % Y value
    %             end
    %             coord2 = coord;
    %             coordD(5,:) = coordD(1,:) ;
    %             coord2(5,:) = coord2(1,:);
    %             subplot(2,2*t2,2*loadcaseIndex-1);
    %             plot(coordD(:,1),coordD(:,2), '-b',coord2(:,1),coord2(:,2), '-g');
    %             axis([-0.3 1.3 -0.3 1.3])
    %             axis square
    %         end
    %     end
    
    % ------------------------------------------------------
    % Generate the Design Vars for the Meso Optimization
    % ------------------------------------------------------
    [DVmeso, configMeso] = GenerateDesignVarsForMesoProblem(config,macroElementProperties.elementNumber ,macroElementProperties);
    
    % ----------------------------------
    % Actual Meso Design.
    % ----------------------------------
    [DVmeso,macroElementProperties,configMeso.totalVolume]= MesoStructureDesignV2(matProp,configMeso,DVmeso,macroElementProperties,[]);
    
    if(configMeso.multiscaleMethodCompare~=1)
        if(config.strainAndTargetTest~=1)
            macroElementProperties.D_subSys  % Show the new D found by the design.
            %   Diff_Sys_Sub =  (macroElementProperties.D_subSys- macroElementProperties.D_sys);
            %     Diff_Sys_Sub
            macroElementProperties.D_sys
            %     determinDiff = det(Diff_Sys_Sub);
            %     determinDiff
            SysDividedbySubSysstem = macroElementProperties.D_sys./macroElementProperties.D_subSys;
            SysDividedbySubSysstem

            %      systemDiff =  macroElementProperties.D_sys-macroElementProperties.D_subSys;
            %     systemDiff
        else
            macroElementProperties.D_subSys 
        end
    end
    
    newDesign = 1;
    if(plottingMesoDesign ==1)
        p = plotResults;
        figure(1)
        if(configMeso.strainAndTargetTest~=1)
            if(configMeso.coordinateMesoBoundaries==1 && configMeso.macro_meso_iteration>1)
                subplot(1,2,1);
                p.PlotArrayGeneric(DVmeso.mesoStructNTCmask,'NTC  sensitivity mask');
            end
            subplot(1,2,2);
            outname = sprintf('meso structure for macro element %i density %f',e, configMeso.totalVolume);
            p.PlotArrayGeneric(DVmeso.x,outname);
               caxis([0 1]);
            %                 subplot(2,2,3);
            %                 outname = sprintf('meso structure sensitivity %i density %f',e, configMeso.v1);
            %                 p.PlotArrayGeneric(DVmeso.temp1,outname);
        else
            pStrain=macroElementProperties.psuedoStrain;
            outname = sprintf('strainAndTargetTest %i density %f',e, configMeso.v1);
            outname2 = sprintf('pseudo strain %f %f %f ',pStrain(1),pStrain(2),pStrain(3));
            outname = [outname outname2];
            p.PlotArrayGeneric(DVmeso.x,outname);
               caxis([0 1]);ng
        end
        drawnow
        nameGraph = sprintf('./out%i/elementpictureIteration%i_element%i.png',configMeso.iterationNum,configMeso.macro_meso_iteration, e);
        print(nameGraph,'-dpng')
    end
else
    %D_homog_flat = zeros(1,9);
    newDesign = 0; % false
    DVmeso=[];
    configMeso=config;
    configMeso.totalVolume=0.1; % just set it to something.
end

SaveMesoUnitCellDesignToCSV(DVmeso,macroElementProperties,configMeso,newDesign);
