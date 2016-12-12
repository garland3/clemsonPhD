function combinedTopologyOptimization(useInputArgs, w1text, macro_meso_iteration,mode, singleMeso_elementNumber)
% useInputArgs = 1 if to use the input args. This expects that we are having a external program call the matlab code
% w1 weight1 for weighted objective.
% macro_meso_iteration, the complete macro-meso iteration number, expect as a string
% mode, if using input args. Then, specify the mode number as a string
% iterationNum, used to know where to output files.
% if mode = 11, then specify the macro element we are designing for. 
close all

% --------------------------------------
%  Settings
% --------------------------------------------
settings = Configuration;

settings.mode =5;
% 1 = topology only,
% 2 = material optimization only.
% 0 = 
% 3 = both mateiral vol fraction and topology
% 4 = meso testing only.
% 5 = Run complete macro optimization, save results to .csv files.
% 6.= Testing reading .csv files and designing each meso structure
% 7= testing recombinging the meso structure designs.
% 10 = macro and meso togehter.
% 11 = single meso design, where "singleMeso_elementNumber" and
% "macro_meso_iteration" specify which one
% 12 Run after meso macro iterations and combines the

macroDesignMode = 0;
if(settings.mode==1 || settings.mode==2 || settings.mode==3 || settings.mode==5 ||  settings.mode==10)
    macroDesignMode=1;
end

settings.parallel = 1;
settings.mesoAddAdjcentCellBoundaries =0;

% material properties Object
matProp = MaterialProperties;

RunPalmetto = 1;
if(RunPalmetto==0)
    % ------------
    % Normal running case
    % -------------------
    settings.macro_meso_iteration = str2double(macro_meso_iteration);
    settings.nelx = 20;
    settings.nely = 20;
    settings.numXElmPerDV=1;
    settings.numYElmPerDV=1;
    
    settings.nelxMeso = 25; %35;
    settings.nelyMeso =25; %35;
    settings.w1 = 1; % do not set to zero, instead set to 0.0001. Else we will get NA for temp2
    settings.iterationNum = 0;
    settings.doSaveDesignVarsToCSVFile = 0;
    settings.doPlotFinal = 1;
    %  settings.terminationCriteria =0.1; % 10%
    settings.terminationCriteria =0.001; % 3%
    settings.numWorkerProcess = 3;
else
    % ------------
    % Palmetto running case
    % -------------------
    settings.nelx = 35;
    settings.nely = 35;
    settings.nelxMeso = 40; %35;
    settings.nelyMeso =40; %35;
    settings.terminationAverageCount = 10;
    settings.terminationCriteria =0.001; % 0.0%
    settings.maxFEACalls = 150;
    settings.maxMasterLoops = 150;
    
end

if(str2double(useInputArgs) ==1)
    % parse potential input arguments.
    settings.macro_meso_iteration=str2double(macro_meso_iteration);
    settings.mode=str2double(mode);
    settings.w1 =str2double(w1text);
    settings= settings.UpdateVolTargetsAndObjectiveWeights();    
    
    % -------------------------
    % Single meso design for MPI parallelization
    % -------------------------
    if ( settings.mode == 11)
        ne = settings.nelx*settings.nely; % number of elements
        elocal =str2double(singleMeso_elementNumber);
        MesoDesignWrapper(settings,elocal, ne,matProp);
        exit
    end
end

settings= settings.UpdateVolTargetsAndObjectiveWeights();
settings


% ---------------------------------
% Mode 12,
% Plot everything together for a particular macro-Meso_iteration
% ---------------------------------
if(settings.mode ==12)
    settings.macro_meso_iteration=5;
    NumMacroMesoIteration= settings.macro_meso_iteration;
    p=plotResults;
    p.PlotEverythingTogether(NumMacroMesoIteration);
end

% ---------------------------------
% Initialization of varriables
% ---------------------------------
designVars = DesignVars(settings);
fractionCurrent_V1Local =0;

if ( macroDesignMode==1)
    % ------------------------------
    % Normal case, 1 design var per element
    % ------------------------------
    designVars.x(1:settings.nely,1:settings.nelx) =settings.v1+settings.v2; % artificial density of the elements
    designVars.w(1:settings.nely,1:settings.nelx)  =settings.v1; % actual volume fraction composition of each element
    
    designVars.sensitivityElastic(1:settings.nely,1:settings.nelx) = 0;
    designVars.sensitivityHeat(1:settings.nely,1:settings.nelx) = 0;
    
    %designVars.complianceSensitivity(1:settings.nely,1:settings.nelx) = 0;
    if (settings.doPlotStress == 1)
        designVars.totalStress(1:settings.nely,1:settings.nelx) = 0;
    end
%     designVars.g1elastic(1:settings.nely,1:settings.nelx) = 0;
%     designVars.g1heat(1:settings.nely,1:settings.nelx) = 0;
    % end
    
    designVars = designVars.CalcIENmatrix(settings);
    designVars =  designVars.CalcNodeLocation(settings);
    designVars = designVars.PreCalculateXYmapToNodeNumber(settings);
    designVars = designVars.CalcElementXYposition(settings);
    
    matProp=  matProp.ReadConstitutiveMatrixesFromFiles(settings);
    
    if (settings.mode == 10 ||settings.mode == 5)
        designVars= ReadXMacroFromCSV( settings,designVars);
    end    
end
masterloop = 0; FEACalls = 0; 



onlyplotfinal =0;
%% ---------------------------------------------------
% Macro Design
% ---------------------------------------------------
% START ITERATION
if(macroDesignMode==1)
    while  masterloop<=settings.maxMasterLoops && FEACalls<=settings.maxFEACalls
        masterloop = masterloop + 1;        
        
        % --------------------------------
        % Run FEA, calculate sensitivities
        % --------------------------------
        designVars = designVars.RunFEAs(settings, matProp, masterloop);
        [vol1Fraction, vol2Fraction] =  designVars.CalculateVolumeFractions(settings);
        
        FEACalls = FEACalls+1;
        
        % --------------------------------
        % Topology Optimization
        % --------------------------------
        if ( settings.mode == 1 || settings.mode == 3 || settings.mode == 10 || settings.mode ==5)      
            designVars = designVars.CalculateTopologySensitivity(settings, matProp, masterloop);
            % normalize the sensitivies  by dividing by their max values.
            if (settings.w1 ~= 1) % if we are using the heat objective
                temp1Max =-1* min(min(designVars.sensitivityElastic));
                designVars.sensitivityElastic = designVars.sensitivityElastic/temp1Max;
                temp2Max = -1* min(min(designVars.sensitivityHeat));
                designVars.sensitivityHeat = designVars.sensitivityHeat/temp2Max;
                
                designVars.dc = settings.w1*designVars.sensitivityElastic+settings.w2*designVars.sensitivityHeat; % add the two sensitivies together using their weights
            else
                designVars.dc = settings.w1*designVars.sensitivityElastic;
            end
            % FILTERING OF SENSITIVITIES            
            [designVars.dc]   = check( settings.nelx, settings.nely,settings.rmin,designVars.x,designVars.dc);
            % DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
            [designVars.x]    = OC( settings.nelx, settings.nely,designVars.x,settings.totalVolume,designVars.dc, designVars, settings);
            % PRINT RESULTS          
            disp([' FEA calls.: ' sprintf('%4i',FEACalls) ' Obj.: ' sprintf('%10.4f',designVars.c) ...
                ' Vol. 1: ' sprintf('%6.3f', vol1Fraction) ...
                ' Vol. 2: ' sprintf('%6.3f', vol2Fraction) ...
                ' Lambda.: ' sprintf('%6.3f',designVars.lambda1  )])
            
            densitySum = sum(sum(designVars.x));
            designVars.storeOptimizationVar = [designVars.storeOptimizationVar;designVars.c, designVars.cCompliance, designVars.cHeat,vol1Fraction,vol2Fraction,fractionCurrent_V1Local,densitySum];
            
            if(onlyplotfinal~=1)
                p = plotResults;
                p.plotTopAndFraction(designVars,  settings, matProp, FEACalls); % plot the results.
            end
           
                % TEST FOR TERMINATION      
            if(TestForTermaination(designVars, settings) ==1)
                disp('break in topology');
                break;
            end
        end % END TOPOLOGY OPTIMIZATION CODE
        
       % --------------------------------
       % Volume fraction optimization
       % --------------------------------
        if ( settings.mode ==2 || settings.mode ==3 || settings.mode == 10 || settings.mode ==5)
            
            designVars = designVars.CalculateMaterialGradientSensitivity(settings, matProp, masterloop);
            
            [vol1Fraction, vol2Fraction] =  designVars.CalculateVolumeFractions(settings);
            
            totalVolLocal = vol1Fraction+ vol2Fraction;
            fractionCurrent_V1Local = vol1Fraction/totalVolLocal;
            targetFraction_v1 = settings.v1/(settings.v1+settings.v2);
            
            % Normalize the sensitives.
            if (settings.w1 ~= 1) % if we are using the heat objective
                temp1Max = max(max(abs(designVars.sensitivityElastic)));
                designVars.sensitivityElastic = designVars.sensitivityElastic/temp1Max;
                temp2Max = max(max(abs(designVars.sensitivityHeat)));
                designVars.sensitivityHeat = designVars.sensitivityHeat/temp2Max;
                
                g1 = settings.w1*designVars.sensitivityElastic+settings.w2*designVars.sensitivityHeat; % Calculate the weighted volume fraction change sensitivity.
            else
                g1 = settings.w1*designVars.sensitivityElastic;
            end
            
            % Filter the g1 sensitivies
            [g1]   = check( settings.nelx, settings.nely,settings.rmin,designVars.x,g1);
            G1 = g1 - designVars.lambda1 +1/(designVars.mu1)*( targetFraction_v1-fractionCurrent_V1Local); % add in the lagrangian
            designVars.w = designVars.w+settings.timestep*G1; % update the volume fraction.            
            designVars.w = max(min( designVars.w,1),0);    % Don't allow the    vol fraction to go above 1 or below 0
            designVars.lambda1 =  designVars.lambda1 -1/(designVars.mu1)*(targetFraction_v1-fractionCurrent_V1Local)*settings.volFractionDamping;
                        
            % PRINT RESULTS           
            densitySum = sum(sum(designVars.x));
            designVars.storeOptimizationVar = [designVars.storeOptimizationVar;designVars.c, designVars.cCompliance, designVars.cHeat,vol1Fraction,vol2Fraction,fractionCurrent_V1Local,densitySum];
            
            if(onlyplotfinal~=1)
                p = plotResults;
                p.plotTopAndFraction(designVars, settings, matProp,FEACalls ); % plot the results.
            end
                        
            disp([' FEA calls.: ' sprintf('%4i',FEACalls) ' Obj.: ' sprintf('%10.4f',designVars.c) ...
                ' Vol. 1: ' sprintf('%6.3f', vol1Fraction) ...
                ' Vol. 2: ' sprintf(    '%6.3f', vol2Fraction) ...
                ' Lambda.: ' sprintf('%6.3f',designVars.lambda1  )])           
           
             % TEST FOR TERMINATION      
            if( TestForTermaination(designVars, settings) ==1)
               disp('break in vol fraction');
                break;
            end
        end % END VOLUME FRACTION OPTIMIZATION CODE
        
        
    end
end

% -----------------------------------------------------------
% PLOT THE FINAL MACRO DESIGN WITH THE STRAINED FEA GRID FOR EACH LOAD
%
% SAVE THE STORED OPTIMIZATION VARRIABLE DATA TO A CSV FILE
% ----------------------------------------------------------
if ( settings.mode ==2 || settings.mode ==3 || settings.mode == 10 || settings.mode ==5)
    folderNum = settings.iterationNum;
    if(1==1)
        [~, t2] = size(settings.loadingCase);        
        for loadcaseIndex = 1:t2
            %   loadcase = settings.loadingCase(loadcaseIndex);
            p = plotResults;
            p.plotTopAndFraction(designVars, settings, matProp,FEACalls ); % plot the results.
            hold on
            p.plotStrainField(settings,designVars,folderNum,loadcaseIndex)
            nameGraph = sprintf('./gradTopOptimization%fwithmesh%i_load%i.png', settings.w1,settings.macro_meso_iteration,loadcaseIndex);
            print(nameGraph,'-dpng');
            hi=  figure(1);
            cla(hi);
            hold off
        end
    end    
    
    outname = sprintf('./out%i/storeOptimizationVarMacroLoop%i.csv',folderNum,settings.macro_meso_iteration);
    csvwrite(outname,designVars.storeOptimizationVar);    
end

% ---------------------------------------------
%
%         MESO DESIGN TESTING, MODE = 4
%         debugging meso design only
%
% ---------------------------------------------
% For testing only
if(settings.mode ==4) % meso-structure design
    TestMesoDesign(designVars,settings,matProp);
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
    SaveMacroProblemStateToCSV(settings,designVars,matProp);
end


%% ---------------------------------------------
%         Meso Design
%
%         LOOP OVER MACRO ELEMENTS AND DESIGN MESO STRUCTURE = MODE 6
%
% ---------------------------------------------
if(settings.mode ==6 ||settings.mode ==8 || settings.mode ==10)
    % the design var object is huge, so I want to garbage collect after
    % saving the important data to disk (.csv files).
    clear designVars
    
    
    %     if(settings.doUseMultiElePerDV==1) % if elements per design var.
    %         ne =  settings.numVarsX*settings.numVarsY;
    %     else
    ne = settings.nelx*settings.nely; % number of elements
    %SavedDmatrix = zeros(ne,9);
    %     end
    
    checkedElements = CalculateCheckedElements(ne, settings);
    allelements = 1:ne;
    nonCheckedElements = setdiff(allelements, checkedElements);
    
    if(settings.parallel==1)
        % Set up parallel computing.
        myCluster = parcluster('local');
        myCluster.NumWorkers = settings.numWorkerProcess;
        saveProfile(myCluster);
        myCluster
        
        poolobj = gcp('nocreate'); % If no pool,create new one.
        if isempty(poolobj)
            parpool('local', settings.numWorkerProcess)
            poolsize = settings.numWorkerProcess;
        else
            poolsize = poolobj.NumWorkers;
        end
        poolsize
        
        % --------------------------------------------------
        % loop over the macro elements and design a meso structure,
        %  parallel using parfor
        % --------------------------------------------------
        parfor_progress(ne);
        
        [~,numElementsInChecked] = size(checkedElements);
        parfor  e = 1:numElementsInChecked
            %             checkedElements = CalculateCheckedElements(ne, settings);
            elocal = checkedElements(e);
            settingscopy = settings; % need for parfor loop.
            settingscopy.useAjacentLocal = 0;
            MesoDesignWrapper(settingscopy,elocal, ne,matProp);
            parfor_progress;
        end
        [~, NumElementsInnonCheckedElements ]= size(nonCheckedElements);
        parfor  e = 1:NumElementsInnonCheckedElements
            elocal = nonCheckedElements(e);
            settingscopy = settings; % need for parfor loop.
            settingscopy.useAjacentLocal = 1;
            MesoDesignWrapper(settingscopy,elocal, ne,matProp);
            parfor_progress;
        end
        parfor_progress(0);
        
    else
        % --------------------------------------------------
        % loop over the macro elements and design a meso structure,
        % no parallel
        % --------------------------------------------------
        
        if(settings.singleMesoDesign~= 1)
            settings.mesoplotfrequency = 50;
            for  e = checkedElements
                settingscopy = settings; % need for parfor loop.
                settingscopy.useAjacentLocal = 0;
                MesoDesignWrapper(settingscopy,e, ne,matProp);
            end
            
            settings.mesoplotfrequency = 50;
            for  e = nonCheckedElements
                settingscopy = settings; % need for parfor loop.
                settingscopy.useAjacentLocal = 1;
                MesoDesignWrapper(settingscopy,e, ne,matProp);
            end
        end
        
        if(settings.singleMesoDesign == 1)
            SingleMesoStuctureWrappe(settings, ne,matProp);
        end
        
    end % end parallel
    clear designVarsMeso
end

% -------------------------------------
% Generate macro-meso complete structure.
% 7 is for testing the recombining of the meso structures.
%
% Loop over the elements and get the design fields, and make one
% huge array showing the actual shape of the structure, tile the
% -------------------------------------
if(settings.mode ==7||settings.mode ==8  || settings.mode ==10 )
    GenerateCompleteStructureV2Improved(settings)
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




function []= GenerateCompleteStructureV2Improved(settings)

postProcess = 1;
close all
p = plotResults;


temp= settings.mesoAddAdjcentCellBoundaries;
settings.mesoAddAdjcentCellBoundaries=0;
[designVarsMeso, mesoSettings] = GenerateDesignVarsForMesoProblem(settings,1);
settings.mesoAddAdjcentCellBoundaries=temp;

mesoSettings.doUseMultiElePerDV =settings.doUseMultiElePerDV;
numTilesX=settings.numTilesX;
numTilesY = settings.numTilesY;

% Generate huge area
totalX=settings.nelx*mesoSettings.nelx*numTilesX
totalY=settings.nely*mesoSettings.nely*numTilesY

completeStruct = zeros(totalY,totalX);
ne = settings.nelx*settings.nely; % number of elements


%--------------------------------------------
% Get the density field
%--------------------------------------------
macro_meso_iteration = settings.macro_meso_iteration;
macroElementProps = macroElementProp;
% macroElementProps.elementNumber = e;
folderNum = settings.iterationNum;
% GET the saved element to XY position map (needed for x and w vars retrival)
outname = sprintf('./out%i/elementXYposition%i.csv',folderNum,macro_meso_iteration);
elementXYposition=csvread(outname);
% Get the density field
outname = sprintf('./out%i/densityfield%i.csv',folderNum,macro_meso_iteration);
xxx = csvread(outname);




for e = 1:ne
    fprintf('element %i of %i\n',e,ne);
    macroElementProps.elementNumber=e;
    results = elementXYposition(macroElementProps.elementNumber,:);
    macroElementProps.yPosition = results(1);
    macroElementProps.xPosition = results(2);
    macroElementProps.density = xxx(macroElementProps.yPosition,macroElementProps.xPosition );
    
    % Check if void
    if(macroElementProps.density>settings.voidMaterialDensityCutOff)
        x=GetMesoUnitCellDesignFromCSV(settings,e);
        designVarsMeso.x = x;
        % -------------------------------------
        % No, post processing
        % -------------------------------------
        if(postProcess~=1)
            yShift = (macroElementProps.yPosition-1)*mesoSettings.nely*numTilesY+1;
            xShift = (macroElementProps.xPosition-1)*mesoSettings.nelx*numTilesX+1;
            designVarsMeso=TileMesoStructure(mesoSettings, designVarsMeso);
            completeStruct(yShift:(yShift+mesoSettings.nely*numTilesY-1),xShift:(xShift+mesoSettings.nelx*numTilesX-1))=designVarsMeso.xTile;
        else
            % -------------------------------------
            % Yes, With, post processing
            % -------------------------------------
            step = 1;
            completeStruct= TileMesoStructureV2(mesoSettings,settings, designVarsMeso,macroElementProps,xxx,completeStruct,step);
            
        end
        
        
    end
end

if(postProcess==1)
    for e = 1:ne
        fprintf('step 2element %i of %i\n',e,ne);
        macroElementProps.elementNumber=e;
        results = elementXYposition(macroElementProps.elementNumber,:);
        macroElementProps.yPosition = results(1);
        macroElementProps.xPosition = results(2);
        macroElementProps.density = xxx(macroElementProps.yPosition,macroElementProps.xPosition );
        
        % Check if void
        %if(macroElementProps.density>settings.voidMaterialDensityCutOff)
        step = 2;
        completeStruct= TileMesoStructureV2(mesoSettings,settings, designVarsMeso,macroElementProps,xxx,completeStruct,step);
        %end
    end
end

% set the max value to be 1
completeStruct( completeStruct>1)=1;

completeStruct(completeStruct>settings.voidMaterialDensityCutOff)=1;
completeStruct(completeStruct<settings.voidMaterialDensityCutOff)=0;


plotname = sprintf('complete structure %i',settings.macro_meso_iteration);
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
nameGraph = sprintf('./completeStucture%f_macroIteration_%i.png', settings.w1,settings.macro_meso_iteration);
print(nameGraph,'-dpng', '-r1200')
outname = sprintf('./completeStucture%f_macroIteration_%i.csv', settings.w1,settings.macro_meso_iteration);
csvwrite(outname,completeStruct);

[xsize, ysize] = size(completeStruct)
ratioLenghtToHieght=2/10;
height = xsize*ratioLenghtToHieght;

% Write ascii STL from gridded data
generateSTL = 0
if(generateSTL==1)
    %     [X,Y] = deal(1:40);             % Create grid reference
    %     X = 1:totalX;
    %     Y = 1:totalY;
    %     Z = completeStruct*height;                  % Create grid height
    %     nameGraph = sprintf('./stl%f_macroIteration_%i.stl', settings.w1,settings.macro_meso_iteration);
    %     stlwrite(nameGraph,X,Y,Z,'mode','binary')
    %   temp(flipud(tril(temp)==temp)) = 0;
    completeStruct(completeStruct>settings.voidMaterialDensityCutOff)=1;
    completeStruct(completeStruct<settings.voidMaterialDensityCutOff)=0;
end
