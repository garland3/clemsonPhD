classdef Configuration
    
    properties
        % --------------------------------------
        % %% Settings
        % --------------------------------------------
        
        % Optimization Configurations
        %
        % each design var will control the density and volume fraction
        % material of several clustered elements.
        mode =4; % 1 = topology only, 2 = material optimization only. 3 = both,4 = meso only. 5 = meso-structure testing
        doUseMultiElePerDV=0; % 1= true, 0 = false do use multiple elements per (1) design varriable.
        averageMultiElementStrain; % 1= true, 0 = false instead of making 1 large strain field, make sevral and average the sensitivies
        singleMesoDesign; % 1 = yes, 0 = true.
        numXElmPerDV= 2; % Number of elements in the X direction for 1 (per) design varriable.
        numYElmPerDV =2; % Number of elements in the X direction for 1 (per) design varriable.
        macro_meso_iteration = 0; % master loop of the whole macro meso system
        mesoAddAdjcentCellBoundaries = 0; % global property
        useAjacentLocal = 0; % local property that maybe be turned on or off
        testingVerGradMaterail = 0;
        numVarsX;
        numVarsY;
        w1 = 1; % weight elastic for multi-objective, % do not set to zero, instead set to 0.0001.
        w2 = 0;
        
        % For ATC optimization
        addConsistencyConstraints=1;      
        Omega = 0.15;
        
        % Modeling settings
        referenceTemperature = 0; % for thermal expansion, assume that there is not strain when at this temperature.
        addThermalExpansion = 0; % Set to 1 to incorporate thermal expansion
        
        % Exx and Eyy Distribution
        useExxEyy=1;       
        useTargetMesoDensity = 1; % 1 = yes, 0 = no and use target Eavg
        targetExxEyyDensity = 0.5;
        useThetaInSurfaceFit = 0;
        
        % rotation
        useRotation =1;
        minRotation =-pi/2;
        maxRotation = pi/2;
        rotationMoveLimit = pi/10;
        
        
        % number of elements
        nelxMeso = 40;
        nelyMeso = 40;
        nelx = 40; % 40 # of elements in the x direcction, must be a multiple of numXElmPerDV
        nely = 20; % 18 number of elements in the y direction, must be a multiple of numYElmPerDV
        
        % SIMP settings
        penal = 3.0; % penality used for the SIMP method
        rmin = 2; % smoothing radius for sensitivity smoothing.
        voidMaterialDensityCutOff = 0.3; % everything below this density is considered void.
        noNewMesoDesignDensityCutOff = 0.28; % any densities below this will not be redesigned. Having a different value than voidMaterialDensityCutOft helps stabalize the algorithm on the structure edges.
        
        
        % VOLUME Fraction SEttings.
        timestep = 0.2; % time step for the volume fraction update algorithm
        volFractionDamping = 0.2;
        v1 = 0.3; % amount of material 1 to use. default to 10%
        v2 = 0.3; % amount of material 2 to use. default to 30%, reduced so there is less meso structures to compute
        totalVolume; % = v1+v2;
        
        
        % Plotting information
        plotSensitivityWhilerunning = 0;
        mesoplotfrequency=100; % how often to plot the meso level design.
        iterationsPerPlot = 5;
        doPlotVolFractionDesignVar = 0;
        doPlotTopologyDesignVar = 0;
        doPlotHeat = 0;
        doPlotHeatSensitivityTopology = 0;
        doPlotStress = 0;
        doPlotFinal =0;
        doPlotMetrics = 0;
        doSaveDesignVarsToCSVFile = 0; % set to 1 to write plotFinal csv file instead
        doPlotAppliedStrain = 0; % For debugging only
        doPlotOrthDistributionVar=0;
        doPlotExx = 0  ;
        doPlotEyy =  0 ;
        doPlotEyyExxArrows =0;
        doPlotElasticSensitivity =  0  ;
        doPlotRotationValue =0;
        doSysANDSubSysDiffValues = 0;
        
        % Exx ,Eyy , Theta (and Rho) Plot Data
        doPlotCombinedExxEyyAndRotation = 1;
        doIncludeRho =1;
        doIncludeSubSystemValues=1;
        %-----------------
        
        recvid = 0; % record video
        maximizePlots= 0;
        
        % ----------------
        % Computational settings
        % ------------------
        iterationNum=0; %  used for parallel computing.
        maxFEACalls = 50;
        maxMasterLoops = 30;
        maxMesoLoops = 80;
        maxNumPseudoStrainLoop=6
        PseudoStrainEndCriteria = 0.1;
        
        terminationAverageCount = 5; % the average change over this number of iterations must be below the termination criteria
        terminationCriteria = 0.001; % if the normalized average change over  terminationAverageCount of iterations is below this value then termainted. ie. 1% change
        % not much faster.
        useGPU = 0; % set to 1 to try to solve matrix using gpu
        parallel =0; % set to 1 to use parfor while preforming the meso design
        numWorkerProcess = 8;
        useCommandLineArgs = 0;
        RunPalmetto =1;
        gr =  (1+sqrt(5))/2 -1 ; % golden ratio  0.618033988749895
        
        % -----------------
        % Use different mixture rules for effective elastic properteis
        % 1. Simple linear interpolation, Vigot rule of miztures E = w(E1)*(1-w)*E2
        % 2. Hashin–Shtrikam law (average of upper and lower boundary)
        % 3. Reuss -rule, 1/E = w/E1+(1-w)/E2 (not implemented yet)
        % 4. Mori and Tanaka, metal ceramic composite
        % ---------------------
        elasticMaterialInterpMethod =1;
        % -----------------
        % Use different mixture rules for effective Heat properteis
        % 1. Simple linear interpolation, Vigot rule of miztures E = w(E1)*(1-w)*E2
        % 3. Reuss -rule, 1/E = w/E1+(1-w)/E2 (not implemented yet)
        % 4. Kingery's, metal ceramic composite
        % 5. Hashin–Shtrikam law (average of upper and lower boundary)
        % ---------------------
        heatMaterialInterpMethod = 1;
        
        % ---------------------------
        % Loading cases
        % ---------------------------
                    loadingCase = [113]; % left clamped, load, middle right
        %           loadingCase = [111 112 113]; % left clamped
        %            loadingCase = [111 112 ]; % left clamped
        %          loadingCase = [111 120 121]; % up, down, right in top right corrner, left clamp.
        %         loadingCase = [111 120]; % up, down, right in top right corrner, left clamp.
        %              loadingCase = [1];
        
        %  loadingCase = [300 301 302 303 304 305]; % shoe
        %  loadingCase = [400 401 402 403 404 405]; % bridge
        %            loadingCase = [404]; % bridge
        %             loadingCase = [113]; % cantilever
        %                 loadingCase = [111]; % top right, force in Y direction
        
        % --------------
        % Meso tiling info
        %--------------
        numTilesX = 3;
        numTilesY = 3;
        
    end
    
    %% -----------------------
    % Methods related to modifying the configuration
    %------------------------
    methods
        function obj = Configuration()
            % Constructor.
        end
        
        
        
        % ------------------------------------
        % Update some of the computed settings
        % -------------------------------------
        function [obj]= UpdateRunTimeComputedSettings(obj, useInputArgs, w1text, macro_meso_iteration,mode, singleMeso_elementNumber)
            
            
            % --------------------------------------
            % PALMETTO
            %
            % If running on the palmetto, then update for palmetto with
            % special settings.
            % --------------------------------------
            if(obj.RunPalmetto==0)
                % ------------
                % Normal running case
                % -------------------
                obj.macro_meso_iteration = str2double(macro_meso_iteration);
                obj.nelx = 20;
                obj.nely = 20;
              
                obj.nelxMeso = 25; %35;
                obj.nelyMeso =25; %35;
                obj.w1 = 1; % do not set to zero, instead set to 0.0001. Else we will get NA for temp2
                obj.iterationNum = 0;
                obj.doSaveDesignVarsToCSVFile = 0;
             
                obj.numWorkerProcess = 3;
            else
                % ------------
                % Palmetto running case
                % -------------------
                obj.nelx = 30; %39 % 30
                obj.nely = 15; % 21 % 15
                obj.nelxMeso = 20; %35;
                obj.nelyMeso =20; %35;
                obj.terminationAverageCount = 10;
                obj.terminationCriteria =0.001; % 0.0%
                obj.maxFEACalls = 200;
                obj.maxMasterLoops = 60;
                
            end
            
            % ----------------------------------------
            % INPUT ARGS
            %
            % Process input arguments
            % -----------------------------------------
            if(str2double(useInputArgs) ==1)
                % parse potential input arguments.
                obj.macro_meso_iteration=str2double(macro_meso_iteration);
                obj.mode=str2double(mode);
                obj.w1 =str2double(w1text);
                
            end
            
            % --------------------------------------
            % Update the volume fraction targets and total volume.
            % --------------------------------------
            % Do not allow w1 to be zero. Divide by zero messes stuff up.
            if(obj.w1 ==0)
                obj.w1=0.00001;
            end
            obj.w2  = 1- obj.w1; % weight heat transfer
            obj.totalVolume = obj.v1+obj.v2;
            
        
        end
        
        
        
        %         function obj = CalculateDesignVarsPerFEAelement(obj)
        %
        %             if(obj.doUseMultiElePerDV ==1)
        %                 if(mod(obj.nely,obj.numYElmPerDV)~=0 || mod(obj.nelx,obj.numXElmPerDV)~=0 )
        %                     disp('nely and numYElmPerDV or nelx and numXElmPerDV not compatible Exiting MATLAB')
        %                     exit
        %                 end
        %
        %                 obj.numVarsX = obj.nelx/obj.numXElmPerDV;
        %                 obj.numVarsY = obj.nely/obj.numYElmPerDV;
        %
        %             end
        %         end
        
        % -----------------------------
        % Given design var position, get list of elements X,Y that are
        % controlled by the design var. Since there are multiple elements
        % the x and y returned values are an array.
        %
        % Only applicable when multiple elements per design var is true.
        %
        % -----------------------------
        %         function [eleXnums, eleYnums,xNodeNums,yNodeNums, macroXdesignVarindex,macroYdesignVarindex] = GetElementsForDesignVar(obj,designvarNumber)
        %             % 1. Get the x,y position of the design var.
        %             % 2. multiply by elements per design var to get positions.
        %
        %             numVarsinRow =obj.numVarsX;
        %             numVarsinColumn =obj.numVarsY;
        %             ydesignVar =   floor( designvarNumber/numVarsinRow)+1;
        %             xdesignVar = mod(designvarNumber/numVarsinColumn);
        %
        %             macroXdesignVarindex = xdesignVar;
        %             macroYdesignVarindex = ydesignVar;
        %
        %             xStart = xdesignVar*obj.numXElmPerDV-(obj.numXElmPerDV);
        %             yStart = ydesignVar*obj.numYElmPerDV-(obj.numYElmPerDV);
        %
        %             numElementsPerDV = obj.numXElmPerDV*obj.numYElmPerDV;
        %             eleXnums = zeros(numElementsPerDV,1);
        %             eleYnums = zeros(numElementsPerDV,1);
        %
        %             count = 1;
        %             for j = 1:obj.numYElmPerDV
        %                 for i = 1:obj.numXElmPerDV
        %                     eleXnums(count) = xStart+i;
        %                     eleYnums(count) = yStart+j;
        %
        %                     count =count+1;
        %                 end
        %             end
        %
        %
        %
        %         end
    end
end