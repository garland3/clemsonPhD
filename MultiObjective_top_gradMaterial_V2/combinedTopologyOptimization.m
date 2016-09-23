function combinedTopologyOptimization(useInputArgs, w1text, iterationNum,macro_meso_iteration)
% input args, useInputArgs = 1 if to use the input args
% w1 weight1 for weighted objective.
% iterationNum, used to know where to output files.
close all
% h = figure(1);hh = figure(2);
% set(h,'visible','off');
% set(hh,'visible','off');

% --------------------------------------
% %% Settings
% --------------------------------------------
settings = Configuration;
settings.macro_meso_iteration = macro_meso_iteration;
settings.mode =6;
% 1 = topology only,
% 2 = material optimization only.
% 3 = both mateiral vol fraction and topology
% 4 = meso testing only.
% 5 = test saving macro structure info to .csv files.  meso-structure
% 6.= Testing reading .csv files and designing each meso structure
% 7= testing recombinging the meso structure designs.
% 8= modes 6 and 7 combined
% 10 = everything working!!!

settings.elasticMaterialInterpMethod = 2; % Hashin–Shtrikam law (average of upper and lower boundary)
settings.heatMaterialInterpMethod = 5; % Hashin–Shtrikam law (average of upper and lower boundary)

% target volumes of material 1 and 2
settings.v1 = 0.2;
settings.v2 = 0.2;

% This is as the new way, but try to make compatible with old
% Each design var controls several elements. 
settings.doUseMultiElePerDV =1;

% if using input args, then override some configurations.
% if using input args, then running on the cluster, so use high resolution,
% otherwise use low resolution
if(str2num(useInputArgs) ==1)
    settings.w1 = str2num(w1text);
    
    settings.iterationNum = str2num(iterationNum);
    settings.nelx = 80;
    settings.nely = 40;
    
    settings.doPlotVolFractionDesignVar = 0;
    settings.doPlotTopologyDesignVar = 0;
    settings.doPlotHeat = 0;
    settings.doPlotHeatSensitivityTopology = 0;
    settings.doPlotStress = 0;
    settings.doPlotFinal = 0;
    settings.doSaveDesignVarsToCSVFile = 1; % set to 1 to write to csv file instead
    settings.maxFEACalls = 350;
    settings.maxMasterLoops = 500; % make it so, the fea maxes out first.
else
    
    settings.nelx = 15;
    settings.nely = 20;
    settings.w1 = 1; % do not set to zero, instead set to 0.0001. Else we will get NA for temp2
    settings.iterationNum = 0;
    settings.doSaveDesignVarsToCSVFile = 0;
    settings.doPlotFinal = 1;
    settings.terminationCriteria =0.1; % 10%
    
end

settings= settings.UpdateVolTargetsAndObjectiveWeights();
settings
onlyTopChangeOnFirstIteration = 1; % 1 = true, 0 = false;
% material properties Object
matProp = MaterialProperties;

% ---------------------------------
% Initialization of varriables
% ---------------------------------
designVars = DesignVars(settings);
designVars.x(1:settings.nely,1:settings.nelx) = settings.totalVolume; % artificial density of the elements
designVars.w(1:settings.nely,1:settings.nelx)  = 1; % actual volume fraction composition of each element
fractionCurrent_V1Local =1;

designVars.temp1(1:settings.nely,1:settings.nelx) = 0;
designVars.temp2(1:settings.nely,1:settings.nelx) = 0;
%designVars.complianceSensitivity(1:settings.nely,1:settings.nelx) = 0;
 
if (settings.doPlotStress == 1)
    designVars.totalStress(1:settings.nely,1:settings.nelx) = 0;
end

designVars.g1elastic(1:settings.nely,1:settings.nelx) = 0;
designVars.g1heat(1:settings.nely,1:settings.nelx) = 0;

designVars = designVars.CalcIENmatrix(settings);
designVars =  designVars.CalcElementLocation(settings);
designVars = designVars.PreCalculateXYmapToNodeNumber(settings);
designVars = designVars.CalcElementXYposition(settings);

masterloop = 0;
FEACalls = 0;
status=0;

% macro_meso_iteration = 1;

if ( settings.mode == 1 || settings.mode == 3 || settings.mode == 10 || settings.mode ==5)
    matProp=  matProp.ReadConstitutiveMatrixesFromFiles( macro_meso_iteration, settings);
end

if ( settings.mode == 10 && onlyTopChangeOnFirstIteration ==1 )
   designVars= ReadXMacroFromCSV(macro_meso_iteration, settings,designVars);
end

% START ITERATION
while status == 0  && masterloop<=settings.maxMasterLoops && FEACalls<=settings.maxFEACalls
    masterloop = masterloop + 1;
    
    %     if macro_meso_iteration>1
    
    %     end
    
    % --------------------------------
    % Topology Optimization
    % --------------------------------
    if ( settings.mode == 1 || settings.mode == 3 || settings.mode == 10 || settings.mode ==5)
        
        if(macro_meso_iteration==1 || (onlyTopChangeOnFirstIteration == 0))
            for loopTop = 1:1
                designVars = designVars.CalculateSensitivies(settings, matProp, masterloop,macro_meso_iteration);
                [vol1Fraction, vol2Fraction] =  designVars.CalculateVolumeFractions(settings);

                FEACalls = FEACalls+1;
                % normalize the sensitivies  by dividing by their max values.
                temp1Max =-1* min(min(designVars.temp1));
                designVars.temp1 = designVars.temp1/temp1Max;
                temp2Max = -1* min(min(designVars.temp2));
                designVars.temp2 = designVars.temp2/temp2Max;

                designVars.dc = settings.w1*designVars.temp1+settings.w2*designVars.temp2; % add the two sensitivies together using their weights

                % FILTERING OF SENSITIVITIES
                [designVars.dc]   = check(settings.nelx,settings.nely,settings.rmin,designVars.x,designVars.dc);
                % DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
                [designVars.x]    = OC(settings.nelx,settings.nely,designVars.x,settings.totalVolume,designVars.dc, designVars, settings);
                % PRINT RESULTS
                %change = max(max(abs(designVars.x-designVars.xold)));
                disp([' FEA calls.: ' sprintf('%4i',FEACalls) ' Obj.: ' sprintf('%10.4f',designVars.c) ...
                    ' Vol. 1: ' sprintf('%6.3f', vol1Fraction) ...
                    ' Vol. 2: ' sprintf('%6.3f', vol2Fraction) ...
                    ' Lambda.: ' sprintf('%6.3f',designVars.lambda1  )])

                densitySum = sum(sum(designVars.x));
                designVars.storeOptimizationVar = [designVars.storeOptimizationVar;designVars.c, designVars.cCompliance, designVars.cHeat,vol1Fraction,vol2Fraction,fractionCurrent_V1Local,densitySum];

                p = plotResults;
                p.plotTopAndFraction(designVars,  settings, matProp, FEACalls); % plot the results.
                status = TestForTermaination(designVars, settings);
                if(status ==1)
                    m = 'break in topology'
                    break;
                end
            end
        end
    end
    
    % exit the master loop if we termination criteria are true.
    if(status ==1)
        m = 'exiting master (break in topology)'
        break;
    end
    
    % --------------------------------
    % Volume fraction optimization
    % --------------------------------
    if ( settings.mode ==2 || settings.mode ==3 || settings.mode == 10 || settings.mode ==5)
        for loopVolFrac = 1:1
            designVars = designVars.CalculateSensitivies( settings, matProp, masterloop,macro_meso_iteration);
            FEACalls = FEACalls+1;
            
            % for j = 1:5
            [vol1Fraction, vol2Fraction] =  designVars.CalculateVolumeFractions(settings);
            
            totalVolLocal = vol1Fraction+ vol2Fraction;
            fractionCurrent_V1Local = vol1Fraction/totalVolLocal;
            targetFraction_v1 = settings.v1/(settings.v1+settings.v2);
            
            % Normalize the sensitives.
            temp1Max = max(max(abs(designVars.g1elastic)));
            designVars.g1elastic = designVars.g1elastic/temp1Max;
            temp2Max = max(max(abs(designVars.g1heat)));
            designVars.g1heat = designVars.g1heat/temp2Max;
            
            g1 = settings.w1*designVars.g1elastic+settings.w2*designVars.g1heat; % Calculate the weighted volume fraction change sensitivity.
            G1 = g1 - designVars.lambda1 +1/(designVars.mu1)*( targetFraction_v1-fractionCurrent_V1Local); % add in the lagrangian
            designVars.w = designVars.w+settings.timestep*G1; % update the volume fraction.
            
            designVars.w = max(min( designVars.w,1),0);    % Don't allow the    vol fraction to go above 1 or below 0
            designVars.lambda1 =  designVars.lambda1 -1/(designVars.mu1)*(targetFraction_v1-fractionCurrent_V1Local)*settings.volFractionDamping;
            
            
            % PRINT RESULTS
            %change = max(max(abs(designVars.x-designVars.xold)));
            densitySum = sum(sum(designVars.x));
            designVars.storeOptimizationVar = [designVars.storeOptimizationVar;designVars.c, designVars.cCompliance, designVars.cHeat,vol1Fraction,vol2Fraction,fractionCurrent_V1Local,densitySum];
            p = plotResults;
            p.plotTopAndFraction(designVars, settings, matProp,FEACalls ); % plot the results.
            
            
            disp([' FEA calls.: ' sprintf('%4i',FEACalls) ' Obj.: ' sprintf('%10.4f',designVars.c) ...
                ' Vol. 1: ' sprintf('%6.3f', vol1Fraction) ...
                ' Vol. 2: ' sprintf(    '%6.3f', vol2Fraction) ...
                ' Lambda.: ' sprintf('%6.3f',designVars.lambda1  )])
            % obj.storeOptimizationVar = [obj.storeOptimizationVar;obj.c,  obj.cCompliance, obj.cHeat ];
            status = TestForTermaination(designVars, settings);
            if(status ==1)
                m = 'break in vol fraction'
                break;
            end
        end
        
        
    end
end

if ( settings.mode ==2 || settings.mode ==3 || settings.mode == 10 || settings.mode ==5)
     folderNum = settings.iterationNum;
    if(1==1)
        plotStrainField(settings,designVars,folderNum,macro_meso_iteration)
    end
    
   
    outname = sprintf('./out%i/storeOptimizationVarMacroLoop%i.csv',folderNum,macro_meso_iteration);
    csvwrite(outname,designVars.storeOptimizationVar);
    status
end
% ---------------------------------------------
%
%         MESO DESIGN TESTING, MODE = 4
%         debugging meso design only
%
% ---------------------------------------------
% For testing only
if(settings.mode ==4) % meso-structure design
%     designVars = designVars.CalcElementNodeMapmatrixWithPeriodicXandY(settings);
%     designVars =  designVars.CalcNodeLocationMeso(settings);
    
%     settings.maxMasterLoops = 3; % make it small
    
    macroElemProps = macroElementProp;
    macroElemProps.material1Fraction = 0.5;
    macroElemProps.disp = [0         0   0   0   .1   0        0.2        0] ; % make these up for now
%     macroElemProps.disp = [0         0   0   0   .1   0        0.1        0] ; % make these up for now
    %macroElemProps.disp = [   -0.0554   -0.0452   -0.0026   -0.0356         0  0     0 0] ; % make these up for now
    mesoSettings = settings;
    mesoSettings.doPlotAppliedStrain = 1;
    MesoStructureDesign(matProp,mesoSettings,designVars,masterloop,FEACalls,macroElemProps);
end

% ---------------------------------------------
%
%         SAVE MACRO PROBLEM TO CSV FILES = MODE 5
%         the macro problem must actually run, so the macro topology and
%         vol fraction must also have an  || settings.mode ==5
%
% ---------------------------------------------

if(settings.mode ==5  || settings.mode ==10)
    
    % write the displacement field to a .csv
    SaveMacroProblemStateToCSV(settings,designVars,macro_meso_iteration,matProp);
end


% ---------------------------------------------
%
%         LOOP OVER MACRO ELEMENTS AND DESIGN MESO STRUCTURE = MODE 6
%
% ---------------------------------------------
if(settings.mode ==6 ||settings.mode ==8 || settings.mode ==10)
    % the design var object is huge, so I want to garbage collect after
    % saving the important data to disk (.csv files).
    clear designVars
    
    ne = settings.nelx*settings.nely; % number of elements
    %SavedDmatrix = zeros(ne,9);
    coord(:,1) = [0 1 1 0];
    coord(:,2)  = [0 0 1 1];
    
    for e = 1:ne
        e
        macroElementProps = GetMacroElementPropertiesFromCSV(settings,e,macro_meso_iteration);
        
        scalePlot = 10;
        
%         plotting = 1; % this was for debugging

           % Check if void
         if(macroElementProps.density>settings.voidMaterialDensityCutOff)
             
              % Only plot a few
              plotfrequency = 1;
                if(mod(e,plotfrequency) ==0)
                      plottingMesoDesign = 1; % this was for debugging
                       plotting = 1; % this was for debugging
                       settings.doPlotAppliedStrain = 1;
                else
                     plottingMesoDesign = 0; % this was for debugging
                       plotting = 0; % this was for debugging
                         settings.doPlotAppliedStrain = 0;
                end
             
            if(plotting ==1)
                figure(1)
                U2 = macroElementProps.disp*scalePlot;
                coordD = zeros(5,2);
                for temp = 1:4
                    coordD(temp,1) =  coord(temp,1)+ U2(2*temp-1); % X value
                    coordD(temp,2) =  coord(temp,2)+ U2(2*temp); % Y value
                end
                coord2 = coord;
                coordD(5,:) = coordD(1,:) ;
                coord2(5,:) = coord2(1,:);
    %             coord3 = coord2;

                subplot(2,2,1);
                plot(coordD(:,1),coordD(:,2), '-b',coord2(:,1),coord2(:,2), '-g');
                axis([-0.3 1.3 -0.3 1.3])
                axis square
            end
        
            [designVarsMeso, mesoSettings] = GenerateDesignVarsForMesoProblem(settings,e,macro_meso_iteration);
            mesoSettings.doPlotAppliedStrain =   settings.doPlotAppliedStrain ;
            
            % Set the target infill for the meso as the vol fraction of
            % material 1 for now.
            %              designVars.x(y,x)
            %             mesoSettings.v1=macroElementProps.material1Fraction*macroElementProps.density^settings.penal;
            mesoSettings.v1=0.5+(macroElementProps.material1Fraction*macroElementProps.density^settings.penal)/2;
            mesoSettings.v2=0;
            mesoSettings.totalVolume= mesoSettings.v1+0;
            [D_homog,designVarsMeso,macroElementProps]= MesoStructureDesign(matProp,mesoSettings,designVarsMeso,masterloop,FEACalls,macroElementProps);
                        D_homog
%             macroElemProps.strain
            %D_homog_flat = reshape(D_homog,1,9);
            
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
                 nameGraph = sprintf('./out%i/elementpicture%i.png',settings.iterationNum, e);
                   print(nameGraph,'-dpng')
            end
            newDesign = 1; % true
        else
            %D_homog_flat = zeros(1,9);
            newDesign = 0; % false
            designVarsMeso=[];
        end
        
        SaveMesoUnitCellDesignToCSV(designVarsMeso,macroElementProps,settings.iterationNum,macro_meso_iteration,e,newDesign);
        
        % SavedDmatrix(e,:) = D_homog_flat;
        
        % write the density, volume fraction and topology fields to .csv files
        % make a list of elemenents that have material (ie, we don't need to
        % design the void regions)
        % loop over the elements with material.
        % 1. read the displacment field and vol frac, calculate the E_given
        % 2. Run the MesoStructureDesign
        % 3. Save the results of the meso-structure to a .csv file.
    end
    
    % Loop over the elements and get the design fields, and make one
    % huge array showing the actual shape of the structure, tile the
    
end

clear designVarsMeso

% -------------------------------------
% Generate macro-meso complete structure.
% 7 is for testing the recombining of the meso structures.
% -------------------------------------
if(settings.mode ==7||settings.mode ==8  || settings.mode ==10 )
    close all
    p = plotResults;
    [designVarsMeso mesoSettings] = GenerateDesignVarsForMesoProblem(settings,1,1);
    % Generate huge area
    numTilesX=settings.numTilesX;
    numTilesY = settings.numTilesY;
    totalX=settings.nelx*mesoSettings.nelx*numTilesX;
    totalY=settings.nely*mesoSettings.nely*numTilesY;
    
    completeStruct = zeros(totalY,totalX);
    ne = settings.nelx*settings.nely; % number of elements
    for e = 1:ne
        e
        macroElementProps = GetMacroElementPropertiesFromCSV(settings,e,macro_meso_iteration);
        
        % Check if void
        if(macroElementProps.density>settings.voidMaterialDensityCutOff)
            x=GetMesoUnitCellDesignFromCSV(settings.iterationNum,macro_meso_iteration,e);
            
            yShift = (macroElementProps.yPosition-1)*mesoSettings.nely*numTilesY+1;
            xShift = (macroElementProps.xPosition-1)*mesoSettings.nelx*numTilesX+1;
            
            designVarsMeso.x = x;
%             subplot(2,2,1);
%              p.PlotArrayGeneric( x, 'local X structure')
             
             designVarsMeso=TileMesoStructure(mesoSettings, designVarsMeso);
%              subplot(2,2,2);
%              p.PlotArrayGeneric( designVarsMeso.xTile, 'local X structure tile')
            
            completeStruct(yShift:(yShift+mesoSettings.nely*numTilesY-1),xShift:(xShift+mesoSettings.nelx*numTilesX-1))=designVarsMeso.xTile;
        end
    end
   
%           subplot(2,2,3);
    p.PlotArrayGeneric( completeStruct, 'complete structure')
    nameGraph = sprintf('./completeStucture%f_macroIteration_%i.png', settings.w1,macro_meso_iteration);
    print(nameGraph,'-dpng')
    
    
end




% test for termaination of the function.
% normalize the objectives, compare to a moving average. If moving average
% below target, then return 1
% if not terminated, then return 0
function status = TestForTermaination(designVars, settings)
status = 0;
y2 = designVars.storeOptimizationVar(:,2); % Elastic Compliance

t = settings.terminationAverageCount;
if(size(y2)<(t+3))
    return;
end


y3 = designVars.storeOptimizationVar(:,3); % Heat Compliance
y4 = designVars.storeOptimizationVar(:,4); % material 1
y5 = designVars.storeOptimizationVar(:,5); % material 2

avg_y2 = FindAvergeChangeOfLastValue(y2,settings); % elastic objective
avg_y3 = FindAvergeChangeOfLastValue(y3,settings); % heat objective
avg_y4 = FindAvergeChangeOfLastValue(y4,settings); % material 1
avg_y5 = FindAvergeChangeOfLastValue(y5,settings); % material 2




tt = settings.terminationCriteria;
if(        avg_y2<tt ...
        && avg_y3<tt ...
        && avg_y4<tt ...
        && avg_y5<tt)
    status=1;
    % print the vars to screen
    [avg_y2 avg_y3 avg_y4 avg_y5 ]
end


function [averageDiffOfLastValues] = FindAvergeChangeOfLastValue(arrayOfValues, settings)
y2 = arrayOfValues;
t = settings.terminationAverageCount;
diff2 = y2-circshift(y2,1); % subtract the current one from the pervious on
diff2 = diff2(2:end-1);

% normalize and take ABS
% we want to normalize by max(abs) over all differences, so
diff2 = abs(diff2)/(max(abs(diff2)));

diff2 = diff2(end-t:end);
averagey2diff = sum(diff2)/t;
averageDiffOfLastValues = averagey2diff;

% if recvid==1
%          close(vidObj);  %close video




%%%%%%%%%% OPTIMALITY CRITERIA UPDATE Gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [wnew]=OC_gradient(g1,designVar, settings)
l1 = 0; l2 = 100000; move = 0.05;
while (l2-l1 > 1e-4)
    lmid = 0.5*(l2+l1);
    designVar.w = max(0.001,max(designVar.w-move,min(1.,min(designVar.w+move,designVar.w.*sqrt(-g1./lmid)))));
    
    %   desvars = max(VOID, max((x - move), min(SOLID,  min((x + move),(x * (-dfc / lammid)**self.eta)**self.q))))
    
    [volume1, volume2] = designVar.CalculateVolumeFractions(settings);
    currentTotal = volume1+volume2;
    %currentvolume=volume1+volume2;
    
    %if currentvolume - volfrac > 0;
    if currentTotal - settings.totalVolume > 0;
        l1 = lmid;
    else
        l2 = lmid;
    end
    
    wnew =  designVar.w ;
end


