function combinedTopologyOptimization(useInputArgs, w1text, macro_meso_iteration,mode, singleMeso_elementNumber)
% useInputArgs = 1 if to use the input args. This expects that we are having a external program call the matlab code
% w1 weight1 for weighted objective.
% macro_meso_iteration, the complete macro-meso iteration number, expect as a string
% mode, if using input args. Then, specify the mode number as a string
% iterationNum, used to know where to output files.
%

% --------------------------------------
%  config
% --------------------------------------------
config = Configuration;
config.mode =65;
% 0- 99 = macro design modes
% 100- 199 = meso design modes
% 200 - 299 = post processing results modes

% 1 = macro topology only,
% 2 = macro material gradient optimization only.
% 3 = macro orthogonal Distribution (old method)
% 4 = macro rotation
% 5 = macro E_xx and E_yy optimization

% 50 = topology,  material gradient
% 55 = topology, material gradient, ortho distribution, rotation
% 60 = topology,  E_xx and E_yy, and rotation together
% 65 = topology,  E_xx and E_yy


% 100 = Single meso design for MPI parallelization
% 110 = TESTING MESO design methods

% 200, Plot the objectives and constraints over several iteraions
% 201, make an .stl file for an iteration.
% 202, Combine meso and Macro designs into a single plot and csv file.
% 203, Extract the Exx, Eyy, Theta, and density values from the meso D matrixes. 


% 3 = both mateiral vol fraction and topology, ortho, rotation
% 4 = meso testing only.
% 5 = Run complete macro optimization, save results to .csv files.
% 6.= Testing reading .csv files and designing each meso structure
% 7= testing recombinging the meso structure designs.
% 10 = macro and meso togehter.
% 11 = single meso design, where "singleMeso_elementNumber" and
% "macro_meso_iteration" specify which one
% 12


opt = Optimizer;

matProp = MaterialProperties; % material properties Object
config.RunPalmetto = 1;
config= config.UpdateRunTimeComputedSettings( useInputArgs, w1text, macro_meso_iteration,mode, singleMeso_elementNumber)


macroDesignMode = 0;
if(config.mode<100)
    macroDesignMode=1;
end

%% -------------------------
% Single meso design for MPI parallelization
% -------------------------
if ( config.mode == 100)
    ne = config.nelx*config.nely; % number of elements
    elocal =str2double(singleMeso_elementNumber);
    MesoDesignWrapper(config,elocal, ne,matProp);
    return
end

%% ---------------------------------------------
%         MESO DESIGN TESTING, MODE = 110
%         debugging meso design only
% ---------------------------------------------
if(config.mode ==110) % meso-structure design
    DV = DesignVars(config);
    TestMesoDesign(DV,config,matProp);
    return;
end

%% ---------------------------------
% Mode 200,
% Plot the objectives and constraints over several iteraions
% ---------------------------------
if(config.mode ==200)
    % config.macro_meso_iteration=5;
    NumMacroMesoIteration= config.macro_meso_iteration;
    p=plotResults;
    p.PlotObjectiveFunctionAndConstraintsOverSevearlIterations(NumMacroMesoIteration,config);
    return
end

%% ---------------------------------
% Mode 201,
% Make an .stl from a .csv file for a particular iteraion.
% ---------------------------------
if(config.mode ==201)
    NumMacroMesoIteration= config.macro_meso_iteration;
    p=plotResults;
    p.PlotEverythingTogether(NumMacroMesoIteration);
    return
end



%% -------------------------------------
% Generate macro-meso complete structure.
% 202 is for  recombining of the meso structures.
% Loop over the elements and get the design fields, and make one
% huge array showing the actual shape of the structure, tile the
% -------------------------------------
if(config.mode ==202)
    GenerateCompleteStructureV2Improved(config)
    return;
end

%% -------------------------------
% Mode 203
% Extract the Exx, Eyy, Theta, and density values from the meso D matrixes
% -------------------------------
if(config.mode ==203)
    GenerateRhoFunctionOfExxEyy(config)
    return;
end


if(config.recvid==1)
    video = VideoManager;
    [vidObj, framedNumber] = video.InitializeVideo( config);
    F=getframe();
end

masterloop = 0; FEACalls = 0;
onlyplotfinal =0;
%% ---------------------------------------------------
% Macro Design
% ---------------------------------------------------
DV = DesignVars(config);
DV.Exx = DV.Exx*config.v1*matProp.E_material1+config.v2*matProp.E_material2;
DV.Eyy = DV.Exx ;

% --------------------------------
% Run FEA, to get started.
% --------------------------------
DV = DV.RunFEAs(config, matProp, masterloop);
DV =  DV.CalculateVolumeFractions(config,matProp);
FEACalls = FEACalls+1;
% START ITERATION
if(macroDesignMode==1)
    while  masterloop<=config.maxMasterLoops && FEACalls<=config.maxFEACalls
        masterloop = masterloop + 1;
        
        % --------------------------------
        % Topology Optimization
        % --------------------------------
        if ( config.mode == 1 || (49 <config.mode  && config.mode < 100))
            DV= opt.OptimizeTopology(DV, config, matProp,masterloop);
            % --------------------------------
            % Run FEA, again.
            % --------------------------------
            DV = DV.RunFEAs(config, matProp, masterloop);
            FEACalls = FEACalls+1;
            DV= DV.CalculateObjectiveValue(config, matProp, masterloop);
            DV = DV.CalculateVolumeFractions(config, matProp) ;
            DV.storeOptimizationVar = [DV.storeOptimizationVar;DV.c, DV.cCompliance, DV.cHeat,DV.currentVol1Fraction,DV.currentVol2Fraction,sum(sum(DV.x)), DV.targetAverageE, DV.actualAverageE];
            
            ShowOptimizerProgress(DV,1,' topology',FEACalls,config, matProp);
            if(TestForTermaination(DV, config,masterloop) ==1)
                disp('break in topology');
                break;
            end
            
        end % END TOPOLOGY OPTIMIZATION CODE
        
        % --------------------------------
        % Volume fraction optimization
        % --------------------------------
        if ( config.mode == 2 ||  config.mode == 50|| config.mode == 55)
            DV= opt.OptimizeVolumeFraction(DV, config, matProp,masterloop);
            % --------------------------------
            % Run FEA, again.
            % --------------------------------
            %             DV = DV.RunFEAs(config, matProp, masterloop);
            %             FEACalls = FEACalls+1;
            DV= DV.CalculateObjectiveValue(config, matProp, masterloop);
            ShowOptimizerProgress(DV,1,' vol fraction',FEACalls,config, matProp);
            DV.storeOptimizationVar = [DV.storeOptimizationVar;DV.c, DV.cCompliance, DV.cHeat,DV.currentVol1Fraction,DV.currentVol2Fraction,sum(sum(DV.x)), DV.targetAverageE, DV.actualAverageE];
            if( TestForTermaination(DV, config,masterloop) ==1)
                disp('break in vol fraction');
                break;
            end
        end % END VOLUME FRACTION OPTIMIZATION CODE
        
        % --------------------------------
        % Rotation of material Optimization
        % --------------------------------
        if(config.useRotation ==1)
            if ( config.mode==4 || config.mode ==55 || config.mode == 60)
                
                DV = opt.OptimizeRotation(DV, config, matProp,masterloop);
              
                % --------------------------------
                % Run FEA,
                % --------------------------------
                %                 DV = DV.RunFEAs(config, matProp, masterloop);
                %                 FEACalls = FEACalls+1;
                DV= DV.CalculateObjectiveValue(config, matProp, masterloop);
                DV.storeOptimizationVar = [DV.storeOptimizationVar;DV.c, DV.cCompliance, DV.cHeat,DV.currentVol1Fraction,DV.currentVol2Fraction,sum(sum(DV.x)), DV.targetAverageE, DV.actualAverageE];
                ShowOptimizerProgress(DV,1,' rotation',FEACalls,config, matProp);
               
                
            end %END ORTHOGONAL MATERIAL DISTRIBUTION OPTIMZATION
        end
        
        % --------------------------------
        % E_xx and E_yy Optimization
        % --------------------------------
        if ( config.mode==5 || config.mode ==60 || config.mode == 65)
            DV = opt.OptimizeExxEyy(DV, config, matProp,masterloop);
            % --------------------------------
            % Run FEA, calculate sensitivities
            % --------------------------------
            %             DV = DV.RunFEAs(config, matProp, masterloop);
            %             FEACalls = FEACalls+1;
            DV = DV.CalculateVolumeFractions(config, matProp) ;
            DV= DV.CalculateObjectiveValue(config, matProp, masterloop);
            DV.storeOptimizationVar = [DV.storeOptimizationVar;DV.c, DV.cCompliance, DV.cHeat,DV.currentVol1Fraction,DV.currentVol2Fraction,sum(sum(DV.x)), DV.targetAverageE, DV.actualAverageE];
            ShowOptimizerProgress(DV,1,' E_xx and E_yy ',FEACalls,config, matProp);
        end % E_xx and E_yy Optimization
        
        if(config.recvid==1)
            [framedNumber, F]  = video.RecordFrame(config,framedNumber, F,vidObj);
        end
    end % MASTER LOOP FOR MACRO LEVEL
    
    % Flip orientation of Exx and Eyy so that theta is positive
    if(config.useRotation ==1)
            if ( config.mode==4 || config.mode ==55 || config.mode == 60)
                 DV = DV.FlipOrientation(config);
            end
    end
    
end % ENDIF FOR MACRO DESIGN


%% -----------------------------------------------------------
% PLOT THE FINAL MACRO DESIGN WITH THE STRAINED FEA GRID FOR EACH LOAD
%
% SAVE THE STORED OPTIMIZATION VARRIABLE DATA TO A CSV FILE
% ----------------------------------------------------------
if ( config.mode <100)
    folderNum = config.iterationNum;
    if(1==1)
        p = plotResults;
        
        % Plot without the displacement Mesh
         p.plotTopAndFraction(DV, config, matProp,FEACalls ); % plot the results.
         nameGraph = sprintf('./gradTopOptimization%fNOhmesh%i.png', config.w1,config.macro_meso_iteration);
         print(nameGraph,'-dpng');
         hi=  figure(1);
         cla(hi);
          % Plot WITH the displacement Mesh
        [~, t2] = size(config.loadingCase);
        for loadcaseIndex = 1:t2               
            p.plotTopAndFraction(DV, config, matProp,FEACalls ); % plot the results.
            hold on
            p.plotStrainField(config,DV,folderNum,loadcaseIndex)
            nameGraph = sprintf('./gradTopOptimization%fwithmesh%i_load%i.png', config.w1,config.macro_meso_iteration,loadcaseIndex);
            print(nameGraph,'-dpng');
            hi=  figure(1);
            if(config.recvid==1)
                [framedNumber, F]  = video.RecordFrame(config,framedNumber, F,vidObj);
            end
            cla(hi);
            hold off
            
        end
    end
    
    outname = sprintf('./out%i/storeOptimizationVarMacroLoop%i.csv',folderNum,config.macro_meso_iteration);
    csvwrite(outname,DV.storeOptimizationVar);
end
if(config.recvid==1)
    video.CloseVideo( config, F,vidObj)
end

%% ---------------------------------------------
%
%         SAVE MACRO PROBLEM TO CSV FILES = MODE 5
%         the macro problem must actually run, so the macro topology and
%         vol fraction must also have an  || config.mode ==5
%
% ---------------------------------------------

if(config.mode <100)
    % write the displacement field to a .csv
    SaveMacroProblemStateToCSV(config,DV,matProp);
end


%% ---------------------------------------------
%         Meso Design
%
%         LOOP OVER MACRO ELEMENTS AND DESIGN MESO STRUCTURE = MODE 6
%
% ---------------------------------------------
% if(config.mode ==6 ||config.mode ==8 || config.mode ==10)
%     % the design var object is huge, so I want to garbage collect after
%     % saving the important data to disk (.csv files).
%     clear DV
%
%
%     %     if(config.doUseMultiElePerDV==1) % if elements per design var.
%     %         ne =  config.numVarsX*config.numVarsY;
%     %     else
%     ne = config.nelx*config.nely; % number of elements
%     %SavedDmatrix = zeros(ne,9);
%     %     end
%
%     checkedElements = CalculateCheckedElements(ne, config);
%     allelements = 1:ne;
%     nonCheckedElements = setdiff(allelements, checkedElements);
%
%     if(config.parallel==1)
%         % Set up parallel computing.
%         myCluster = parcluster('local');
%         myCluster.NumWorkers = config.numWorkerProcess;
%         saveProfile(myCluster);
%         myCluster
%
%         poolobj = gcp('nocreate'); % If no pool,create new one.
%         if isempty(poolobj)
%             parpool('local', config.numWorkerProcess)
%             poolsize = config.numWorkerProcess;
%         else
%             poolsize = poolobj.NumWorkers;
%         end
%         poolsize
%
%         % --------------------------------------------------
%         % loop over the macro elements and design a meso structure,
%         %  parallel using parfor
%         % --------------------------------------------------
%         parfor_progress(ne);
%
%         [~,numElementsInChecked] = size(checkedElements);
%         parfor  e = 1:numElementsInChecked
%             %             checkedElements = CalculateCheckedElements(ne, config);
%             elocal = checkedElements(e);
%             configcopy = config; % need for parfor loop.
%             configcopy.useAjacentLocal = 0;
%             MesoDesignWrapper(configcopy,elocal, ne,matProp);
%             parfor_progress;
%         end
%         [~, NumElementsInnonCheckedElements ]= size(nonCheckedElements);
%         parfor  e = 1:NumElementsInnonCheckedElements
%             elocal = nonCheckedElements(e);
%             configcopy = config; % need for parfor loop.
%             configcopy.useAjacentLocal = 1;
%             MesoDesignWrapper(configcopy,elocal, ne,matProp);
%             parfor_progress;
%         end
%         parfor_progress(0);
%
%     else
%         % --------------------------------------------------
%         % loop over the macro elements and design a meso structure,
%         % no parallel
%         % --------------------------------------------------
%
%         if(config.singleMesoDesign~= 1)
%             config.mesoplotfrequency = 50;
%             for  e = checkedElements
%                 configcopy = config; % need for parfor loop.
%                 configcopy.useAjacentLocal = 0;
%                 MesoDesignWrapper(configcopy,e, ne,matProp);
%             end
%
%             config.mesoplotfrequency = 50;
%             for  e = nonCheckedElements
%                 configcopy = config; % need for parfor loop.
%                 configcopy.useAjacentLocal = 1;
%                 MesoDesignWrapper(configcopy,e, ne,matProp);
%             end
%         end
%
%         if(config.singleMesoDesign == 1)
%             SingleMesoStuctureWrappe(config, ne,matProp);
%         end
%
%     end % end parallel
%     clear DVMeso
% end








% test for termaination of the function.
% normalize the objectives, compare to a moving average. If moving average
% below target, then return 1
% if not terminated, then return 0
function status = TestForTermaination(DV, config,masterloop)

if(masterloop<3)
    status=0;
    return;
end
status = 0;
y2 = DV.storeOptimizationVar(:,2); % Elastic Compliance

t = config.terminationAverageCount;
if(size(y2)<(t+3))
    return;
end


y3 = DV.storeOptimizationVar(:,3); % Heat Compliance
y4 = DV.storeOptimizationVar(:,4); % material 1
y5 = DV.storeOptimizationVar(:,5); % material 2

avg_y2 = FindAvergeChangeOfLastValue(y2,config); % elastic objective
avg_y3 = FindAvergeChangeOfLastValue(y3,config); % heat objective
avg_y4 = FindAvergeChangeOfLastValue(y4,config); % material 1
avg_y5 = FindAvergeChangeOfLastValue(y5,config); % material 2




tt = config.terminationCriteria;
if(        avg_y2<tt ...
        && avg_y3<tt ...
        && avg_y4<tt ...
        && avg_y5<tt)
    status=1;
    % print the vars to screen
    [avg_y2 avg_y3 avg_y4 avg_y5 ]
end

% -----------------------
% SHOW STATUS OF OPTIMIZER
% Print and Plot information.
% Store data in the storeOptimizationVar
% -----------------------
function ShowOptimizerProgress(DV,doPlot,name,FEACalls,config, matProp)


% PRINT RESULTS
disp([' FEA calls.: ' sprintf('%4i',FEACalls) ' Obj.: ' sprintf('%10.4f',DV.c) ' Vol. 1: ' sprintf('%6.3f', DV.currentVol1Fraction)  ' Vol. 2: ' sprintf('%6.3f', DV.currentVol2Fraction) ...
    ' Target E.: ' sprintf('%4i',DV.targetAverageE)    ' Current E.: ' sprintf('%4i',DV.actualAverageE) name])

if(doPlot==1)
    p = plotResults;
    p.plotTopAndFraction(DV,  config, matProp, FEACalls); % plot the results.
end

function []= GenerateCompleteStructureV2Improved(config)

postProcess = 1;
close all
p = plotResults;


temp = config.mesoAddAdjcentCellBoundaries;
config.mesoAddAdjcentCellBoundaries=0;
macroElementProps = macroElementProp;
macroElementProps.targetDensity=0.5; % make up some value for now. 
[DVMeso, mesoconfig] = GenerateDesignVarsForMesoProblem(config,1,macroElementProps);
config.mesoAddAdjcentCellBoundaries=temp;

mesoconfig.doUseMultiElePerDV =config.doUseMultiElePerDV;
numTilesX=config.numTilesX;
numTilesY = config.numTilesY;

% Generate huge area
totalX=config.nelx*mesoconfig.nelx*numTilesX
totalY=config.nely*mesoconfig.nely*numTilesY

completeStruct = zeros(totalY,totalX);
ne = config.nelx*config.nely; % number of elements


%--------------------------------------------
% Get the density field
%--------------------------------------------
macro_meso_iteration = config.macro_meso_iteration;
%macroElementProps = macroElementProp;
% macroElementProps.elementNumber = e;
folderNum = config.iterationNum;
% GET the saved element to XY position map (needed for x and w vars retrival)
outname = sprintf('./out%i/elementXYposition%i.csv',folderNum,macro_meso_iteration);
elementXYposition=csvread(outname);
% Get the density field
outname = sprintf('./out%i/SIMPdensityfield%i.csv',folderNum,macro_meso_iteration);
xxx = csvread(outname);




for e = 1:ne
    fprintf('element %i of %i\n',e,ne);
    macroElementProps.elementNumber=e;
    results = elementXYposition(macroElementProps.elementNumber,:);
    macroElementProps.yPos = results(1);
    macroElementProps.xPos = results(2);
    macroElementProps.densitySIMP = xxx(macroElementProps.yPos,macroElementProps.xPos );
    
    % Check if void
    if(macroElementProps.densitySIMP>config.voidMaterialDensityCutOff)
        x=GetMesoUnitCellDesignFromCSV(config,e);
        DVMeso.x = x;
        % -------------------------------------
        % No, post processing
        % -------------------------------------
        if(postProcess~=1)
            yShift = (macroElementProps.yPos-1)*mesoconfig.nely*numTilesY+1;
            xShift = (macroElementProps.xPos-1)*mesoconfig.nelx*numTilesX+1;
            DVMeso=TileMesoStructure(mesoconfig, DVMeso);
            completeStruct(yShift:(yShift+mesoconfig.nely*numTilesY-1),xShift:(xShift+mesoconfig.nelx*numTilesX-1))=DVMeso.xTile;
        else
            % -------------------------------------
            % Yes, With, post processing
            % -------------------------------------
            step = 1;
            completeStruct= TileMesoStructureV2(mesoconfig,config, DVMeso,macroElementProps,xxx,completeStruct,step);
            
        end
        
        
    end
end

if(postProcess==1)
    for e = 1:ne
        fprintf('step 2element %i of %i\n',e,ne);
        macroElementProps.elementNumber=e;
        results = elementXYposition(macroElementProps.elementNumber,:);
        macroElementProps.yPos = results(1);
        macroElementProps.xPos = results(2);
        macroElementProps.densitySIMP = xxx(macroElementProps.yPos,macroElementProps.xPos );
        
        % Check if void
        %if(macroElementProps.densitySIMP>config.voidMaterialDensityCutOff)
        step = 2;
        completeStruct= TileMesoStructureV2(mesoconfig,config, DVMeso,macroElementProps,xxx,completeStruct,step);
        %end
    end
end

% set the max value to be 1
completeStruct( completeStruct>1)=1;

completeStruct(completeStruct>config.voidMaterialDensityCutOff)=1;
completeStruct(completeStruct<config.voidMaterialDensityCutOff)=0;


plotname = sprintf('complete structure %i',config.macro_meso_iteration);
p.PlotArrayGeneric( completeStruct, plotname)
rgbSteps = 100;  caxis([0,1]);
map = colormap; % current colormap
middPoint = floor(rgbSteps/4);
map(1:middPoint,:) = [ones(middPoint,1),ones(middPoint,1),ones(middPoint,1)];
for zz =    middPoint:rgbSteps
    map(zz,:) = [0,               1- zz/rgbSteps, 0.5];
end
colormap(map)
%     colorbar
freezeColors
nameGraph = sprintf('./completeStucture%f_macroIteration_%i.png', config.w1,config.macro_meso_iteration);
print(nameGraph,'-dpng', '-r1200')
outname = sprintf('./completeStucture%f_macroIteration_%i.csv', config.w1,config.macro_meso_iteration);
csvwrite(outname,completeStruct);
