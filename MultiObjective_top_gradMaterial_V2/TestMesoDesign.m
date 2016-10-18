function TestMesoDesign(designVars,settings,matProp)

macroElemProps = macroElementProp;
macroElemProps.material1Fraction = 0.5;
scalePlot = 30;

e = 1;
macroElemProps = GetMacroElementPropertiesFromCSV(settings,e);
[designVars, settings] = GenerateDesignVarsForMesoProblem(settings,e);

coord(:,1) = [0 1 1 0];
coord(:,2)  = [0 0 1 1];
    settings.loadingCase = [1];
    
settings.doPlotAppliedStrain = 1;  plottingMesoDesign = 1;  plotting = 1; % this was for debugging % this was for debugging
% plotting= 1;
if(settings.doUseMultiElePerDV==0)
    if(plotting ==1)
        figure(1)
        macroElemProps.disp = [0         0   0   0   .1   0        0.2        0] ; % make these up for now
    
        [~, t2] = size(settings.loadingCase);
        for loadcaseIndex = 1:t2
            % utemp = U(loadcaseIndex,:);
            U2 = macroElemProps.disp(loadcaseIndex,:)*scalePlot;
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
     macroElemProps.xDisplacements = [ 0 0 0 0 0 0  1 1 1]*0.01;
     macroElemProps.yDisplacements = [ 0 0 0 0 0 0  1 1 1]*0.01;
       
       
    [Y,X] = ndgrid(0:settings.numYElmPerDV,0:settings.numXElmPerDV);
    [t1, t2] = size(X);

    
    macroElemProps.mesoXnodelocations=X;
    macroElemProps.mesoYnodelocations=Y;
    
    
    if(plotting ==1)
        figure(1)
        [~, t2] = size(settings.loadingCase);
        for loadcaseIndex = 1:t2
            dx = macroElemProps.xDisplacements(loadcaseIndex,:)*scalePlot;
            dy = macroElemProps.yDisplacements(loadcaseIndex,:)*scalePlot;
            Xlocs = macroElemProps.mesoXnodelocations;
            Ylocs = macroElemProps.mesoYnodelocations;
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


%     macroElemProps.disp = [0         0   0   0   .1   0        0.1        0] ; % make these up for now
%macroElemProps.disp = [   -0.0554   -0.0452   -0.0026   -0.0356         0  0     0 0] ; % make these up for now
mesoSettings = settings;
mesoSettings.doPlotAppliedStrain = 1;
masterloop = 1;
FEACalls = 1;
dcGiven = 0;
mesoSettings.doPlotAppliedStrain = 1;  
MesoStructureDesign(matProp,mesoSettings,designVars,masterloop,FEACalls,macroElemProps,dcGiven);