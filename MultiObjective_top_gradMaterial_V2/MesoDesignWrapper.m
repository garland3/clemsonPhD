function  MesoDesignWrapper(settingscopy,e,ne,matProp)

disp(['Meso Design #: ' sprintf('%4i',e ) ' of ' sprintf('%4i',ne )]);
macroElementPropsParFor = GetMacroElementPropertiesFromCSV(settingscopy,e);
scalePlot = 1;
coord(:,1) = [0 1 1 0];
coord(:,2)  = [0 0 1 1];

% Check if void
%         if(macroElementPropsParFor.density>settingscopy.voidMaterialDensityCutOff)
if(macroElementPropsParFor.density>0.05)
    
    % Only plot a few
    settingscopy.mesoplotfrequency = 200;
    if(mod(e,settingscopy.mesoplotfrequency) ==0)
        settingscopy.doPlotAppliedStrain = 1;  plottingMesoDesign = 1;  plotting = 1; % this was for debugging % this was for debugging
    else
        settingscopy.doPlotAppliedStrain = 0; plottingMesoDesign = 0;    plotting = 0; % this was for debugging
    end
    
    % ------------------------------------------------------
    % Single element per design var.
    % ------------------------------------------------------
    if(settingscopy.doUseMultiElePerDV==0)
        if(plotting ==1)
            figure(1)
            [~, t2] = size(settingscopy.loadingCase);
            for loadcaseIndex = 1:t2
                % utemp = U(loadcaseIndex,:);
                U2 = macroElementPropsParFor.disp(loadcaseIndex,:)*scalePlot;
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
        
    else
        % ------------------------------------------------------
        % Multiple element per design var.
        % ------------------------------------------------------
        
        if(e==-1)
            macroElementPropsParFor.xDisplacements = [ 0 0 0 0 0 0  1 1 1]*0.1;
            macroElementPropsParFor.yDisplacements = [ 0 0 0 0 0 0  1 1 1]*0.1;
        end
        
        if(plotting ==1)
            figure(1)
            [~, t2] = size(settingscopy.loadingCase);
            for loadcaseIndex = 1:t2
                dx = macroElementPropsParFor.xDisplacements(loadcaseIndex,:)*scalePlot;
                dy = macroElementPropsParFor.yDisplacements(loadcaseIndex,:)*scalePlot;
                Xlocs = macroElementPropsParFor.mesoXnodelocations;
                Ylocs = macroElementPropsParFor.mesoYnodelocations;
                Xlocs = reshape(Xlocs',[],1);
                Ylocs = reshape(Ylocs',[],1);
                displacedX= Xlocs +dx';
                displacedY= Ylocs +dy';
                subplot(2,2*t2,2*loadcaseIndex-1);
                plot(Xlocs,Ylocs,'x');
                hold on
                plot(displacedX,displacedY,'o');
                hold off
            end
        end
    end
    
    [designVarsMeso, mesoSettings] = GenerateDesignVarsForMesoProblem(settingscopy,e);
    
    
    % Set the target infill for the meso as the vol fraction of
    %             mesoSettings.v1=0.5+(macroElementPropsParFor.material1Fraction*macroElementPropsParFor.density^settings.penal)/2;
    w = macroElementPropsParFor.material1Fraction;
    x = macroElementPropsParFor.density;
    mesoSettings.v1=matProp.CalculateDensityTargetforMeso(w,x,settingscopy);
    mesoSettings.v2=0;
    mesoSettings.totalVolume= mesoSettings.v1+0;
    
    
    mesoSettings.averageMultiElementStrain= settingscopy.averageMultiElementStrain;
    mesoSettings.doPlotAppliedStrain=settingscopy.doPlotAppliedStrain;
    [D_homog,designVarsMeso,macroElementPropsParFor]= MesoStructureDesignV2(matProp,mesoSettings,designVarsMeso,macroElementPropsParFor,[]);
    D_homog
    
    newDesign = 1;
    
    
    if(plottingMesoDesign ==1)
        p = plotResults;
        figure(1)
        subplot(2,2,2);
        outname = sprintf('meso structure for macro element %i density %f',e, mesoSettings.v1);
        p.PlotArrayGeneric(designVarsMeso.x,outname);
        %                 subplot(2,2,3);
        %                 outname = sprintf('meso structure sensitivity %i density %f',e, mesoSettings.v1);
        %                 p.PlotArrayGeneric(designVarsMeso.temp1,outname);
        drawnow
        nameGraph = sprintf('./out%i/elementpicture%i.png',settingscopy.iterationNum, e);
        print(nameGraph,'-dpng')
    end
else
    %D_homog_flat = zeros(1,9);
    newDesign = 0; % false
    designVarsMeso=[];
end

SaveMesoUnitCellDesignToCSV(designVarsMeso,macroElementPropsParFor,settingscopy.iterationNum,settingscopy.macro_meso_iteration,e,newDesign);

% SavedDmatrix(e,:) = D_homog_flat;

% write the density, volume fraction and topology fields to .csv files
% make a list of elemenents that have material (ie, we don't need to
% design the void regions)
% loop over the elements with material.
% 1. read the displacment field and vol frac, calculate the E_given
% 2. Run the MesoStructureDesign
% 3. Save the results of the meso-structure to a .csv file.