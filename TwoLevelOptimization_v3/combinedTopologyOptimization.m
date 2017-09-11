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
% 6 = anisotropic material


% 50 = topology,  material gradient
% 55 = topology, material gradient, ortho distribution, rotation
% 60 = topology,  E_xx and E_yy, and rotation together
% 65 = topology,  E_xx and E_yy
% 70 = topology and anisotropicMat
% 80 = testing response surface for Exx,Eyy, rho
% 90,290, Calculate objective of macro using D sub system matrixes. 


% 100 = Single meso design for MPI parallelization
% 110 = TESTING MESO design methods
% 111 = Validate Meso (generate Targets)
% 112 = Read the meso design information and compute validation metrics
% 113 = Generate Pseudo Strain and Density Targts Experimentt
% 114 = Read Psuedo strain and density target results. Save .csv file
% 115 = interprete the psuedo strain and density target results as graphs

% 200, Plot the objectives and constraints over several iteraions
% 201, make an .stl file for an iteration.
% 202, Combine meso and Macro designs into a single plot and csv file.
% 203, Extract the Exx, Eyy, Theta, and density values from the meso D matrixes. 
% 204, ANN
% 


opt = Optimizer;

matProp = MaterialProperties; % material properties Object
config= config.UpdateRunTimeComputedSettings( useInputArgs, w1text, macro_meso_iteration,mode, singleMeso_elementNumber);


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

%% ---------------------------------------------
%          Validate Meso (generate Targets), MODE = 111        
% ---------------------------------------------
if(config.mode ==111) 
    DV = DesignVars(config);
    step=1;
    GenerateMesoValidationTargets(DV,config,matProp,step);
    return;
end

%% ---------------------------------------------
%         Read the meso design information, MODE = 112
%        and compute validation metrics 
% ---------------------------------------------
if(config.mode ==112) 
%     DV = DesignVars(config);
    fprintf('compute Meso Design Metrics.\n');
       GenerateRhoFunctionOfExxEyy(config)
%        step=2;
%         DV = DesignVars(config);
%         GenerateMesoValidationTargets(DV,config,matProp,step);
%    ComputeMesoDesignMetrics(DV,config,matProp);
    return;
end


%% ---------------------------------------------
%        Generate Pseudo Strain and Density Targts for ANN map test
% ---------------------------------------------
if(config.mode ==113)     
    fprintf('GeneratePsuedoStrainsAndDensityTargets\n');
     DV = DesignVars(config);
     step=1;
    GeneratePsuedoStrainsAndDensityTargets(DV,config,matProp,step);
    return; 
end

%% ---------------------------------------------
%       Read the  Generated Pseudo Strain and Density Targts for ANN map test
% ---------------------------------------------
if(config.mode ==114 || config.mode ==115)     
   
     DV = DesignVars(config);
     if(config.mode ==114)
          fprintf('GeneratePsuedoStrainsAndDensityTargets STEP 2\n');
        step=2;
         GeneratePsuedoStrainsAndDensityTargets(DV,config,matProp,step);
     elseif(config.mode ==115)
          fprintf('GeneratePsuedoStrainsAndDensityTargets STEP 3\n');
          step=3;
         GeneratePsuedoStrainsAndDensityTargets(DV,config,matProp,step);
         
     end
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
    p.PlotObjectiveFunctionAndConstraintsOverSevearlIterations(NumMacroMesoIteration,config,matProp);
    return
end

%% ---------------------------------
% Mode 201,
% Make an .stl from a .csv file for a particular iteraion.
% ---------------------------------
if(config.mode ==201)
    NumMacroMesoIteration= config.macro_meso_iteration;
%     p=plotResults;
%     p.PlotEverythingTogether(NumMacroMesoIteration);
outname = sprintf('./completeStucture%f_macroIteration_%i.csv', config.w1,config.macro_meso_iteration);
 ConvertCSVToSTL(outname)
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

%% -------------------------------
% Mode 204
% Calculate objective of macro using D sub system matrixes. 
% -------------------------------
if(config.mode ==204)
    annTest(1);
    return;
end

if(config.recvid==1)
    video = VideoManager;
    [vidObj, framedNumber] = video.InitializeVideo( config,'macroOptVideoTopology.avi');
    F=getframe();
end

masterloop = 0; FEACalls = 0;

%% ---------------------------------------------------
% Macro Design
% ---------------------------------------------------
DV = DesignVars(config);

if(config.mode==90)
    matProp= matProp.ReadConstitutiveMatrixesFromFiles(config);
end

% if iteration 2 or higher, then get saved problems state and calcualte
% penalty function valeus. 
 if (49<config.mode && config.mode <100  )
    DV=DV.GetMacroStateVarsFromCSVFiles( config,matProp);
    DV=DV.UpdatePenaltyAndLagrangianValues( config,matProp);
 end

 if ( config.mode == 1)
     DV.Exx = ones(size(DV.Exx))*matProp.E_material1;
     DV.Eyy = DV.Exx ;
 end
 
 %  opt=opt.GenerateInterpolateANN(DV.ResponseSurfaceCoefficents,config,matProp);
 if(config.multiscaleMethodCompare==1)
     DV=DV.GetMacroStateVarsFromCSVFiles( config);
    matProp= matProp.ReadConstitutiveMatrixesFromFiles(config);
     % USe in FEA???
 end

% --------------------------------
% Run FEA, to get started.
% --------------------------------
DV = DV.RunFEAs(config, matProp, masterloop);
DV =  DV.CalculateVolumeFractions(config,matProp);
FEACalls = FEACalls+1;
% START ITERATION
if(macroDesignMode==1 &&  config.mode ~= 90)
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
            if(config.multiscaleMethodCompare~=1)
                DV = DV.RunFEAs(config, matProp, masterloop);               
            end
           
            DV= DV.CalculateObjectiveValue(config, matProp, masterloop,opt);
            DV = DV.CalculateVolumeFractions(config, matProp) ;
            DV = DV.AddDataToStoreOptimizationVarArray(config);
            ShowOptimizerProgress(DV,1,' topology',FEACalls,config, matProp);
             FEACalls = FEACalls+1;
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
            DV= DV.CalculateObjectiveValue(config, matProp, masterloop,opt);
             DV = DV.AddDataToStoreOptimizationVarArray(config);
            ShowOptimizerProgress(DV,1,' vol fraction',FEACalls,config, matProp);
           
            if( TestForTermaination(DV, config,masterloop) ==1)
                disp('break in vol fraction');
                break;
            end
        end % END VOLUME FRACTION OPTIMIZATION CODE
        
        % --------------------------------
        % anisotropic material optimization
        % --------------------------------
        if ( config.mode == 6 ||  config.mode == 70)
            DV= opt.OptimizeAnisotropicMaterial(DV, config, matProp,masterloop);
            % --------------------------------
            % Run FEA, again.
            % --------------------------------
            DV = DV.RunFEAs(config, matProp, masterloop);
            FEACalls = FEACalls+1;
            DV= DV.CalculateObjectiveValue(config, matProp, masterloop,opt);
            ShowOptimizerProgress(DV,1,' anisotropic',FEACalls,config, matProp);
            DV = DV.AddDataToStoreOptimizationVarArray(config);
            if( TestForTermaination(DV, config,masterloop) ==1)
                disp('break in anisotropic');
                break;
            end
        end % END anisotropic OPTIMIZATION CODE
        
        % --------------------------------
        % Rotation of material Optimization
        % --------------------------------
        if(config.useRotation ==1)
            if(config.simplifiedOrth_noRotation==0)
            if ( config.mode==4 || config.mode ==55 || config.mode == 60)
                
%                 if(mod(masterloop,5)==1 || config.macro_meso_iteration==1)
                    
                    DV = opt.OptimizeRotation(DV, config, matProp,masterloop);
                    
                    % --------------------------------
                    % Run FEA,
                    % --------------------------------
                    DV = DV.RunFEAs(config, matProp, masterloop);
                 
                    DV= DV.CalculateObjectiveValue(config, matProp, masterloop,opt);
                    DV = DV.AddDataToStoreOptimizationVarArray(config);
                    ShowOptimizerProgress(DV,1,' rotation',FEACalls,config, matProp);
                       FEACalls = FEACalls+1;
                    
%                 end
            end %END ORTHOGONAL MATERIAL DISTRIBUTION OPTIMZATION
            end
        end
        
        % --------------------------------
        % E_xx and E_yy Optimization
        % --------------------------------
        if ( config.mode==5 || config.mode ==60 || config.mode == 65)
            DV = opt.OptimizeExxEyy(DV, config, matProp,masterloop);
            % --------------------------------
            % Run FEA, calculate sensitivities
            % --------------------------------
                        DV = DV.RunFEAs(config, matProp, masterloop);
                      
            DV = DV.CalculateVolumeFractions(config, matProp) ;
            DV= DV.CalculateObjectiveValue(config, matProp, masterloop,opt);
             DV = DV.AddDataToStoreOptimizationVarArray(config);
            ShowOptimizerProgress(DV,1,' E_xx and E_yy ',FEACalls,config, matProp);
           FEACalls = FEACalls+1;
        end % E_xx and E_yy Optimization
        
          % Flip orientation of Exx and Eyy so that theta is positive
        if(config.useRotation ==1)
            if ( config.mode==4 || config.mode ==55 || config.mode == 60)
                DV = DV.FlipOrientation(config);
%                 ShowOptimizerProgress(DV,1,' Flipped theta ',FEACalls,config, matProp);
%                 nameGraph = sprintf('./FinalWithFlip%f__%i.png', config.w1,config.macro_meso_iteration);
%                 print(nameGraph,'-dpng');
            end
        end
              
        
        if(config.recvid==1)
            [framedNumber, F]  = video.RecordFrame(config,framedNumber, F,vidObj);
        end
   
     
      
     end % MASTER LOOP FOR MACRO LEVEL
   %  nameGraph = sprintf('./gradTopOptimizationPReFlip%fNOhmesh%i.png', config.w1,config.macro_meso_iteration);
    % print(nameGraph,'-dpng');

    
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
    
    
end
if(config.recvid==1)
    video.CloseVideo( config, F,vidObj)
end
if(config.mode ==90)
    DV= DV.CalculateObjectiveValue(config, matProp, masterloop,opt);
    DV = DV.CalculateVolumeFractions(config, matProp) ;
%     DV = AddDataToStoreOptimizationVarArray(DV,config);
    %        ShowOptimizerProgress(DV,1,' MODE 90 ',FEACalls,config, matProp);
    outname = sprintf('./mode90/mode90_objective%i.csv',config.macro_meso_iteration);
    csvwrite(outname,DV.c);
     outname = sprintf('./out%i/storeOptimizationVarMacroLoop%i.csv',folderNum,config.macro_meso_iteration);
    storeOptimizationVar=csvread(outname);
    elasticObjectiveMacroModel = storeOptimizationVar(end,2); % Elastic Compliance
    fprintf('Mode 90, final objective value uisng D_sub values %f, compared to %f using macro model\n',DV.c,elasticObjectiveMacroModel);
    plotMode90Data(config.macro_meso_iteration);
end

%% ---------------------------------------------
%
%         SAVE MACRO PROBLEM TO CSV FILES = MODE 5
%         the macro problem must actually run, so the macro topology and
%         vol fraction must also have an  || config.mode ==5
%
% ---------------------------------------------

if(config.mode <100 && config.mode~=90)
    outname = sprintf('./out%i/storeOptimizationVarMacroLoop%i.csv',folderNum,config.macro_meso_iteration);
    csvwrite(outname,DV.storeOptimizationVar);
    % write the displacement field to a .csv
    SaveMacroProblemStateToCSV(config,DV,matProp);
end





