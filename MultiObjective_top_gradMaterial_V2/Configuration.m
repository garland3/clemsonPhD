classdef Configuration
    
    properties
        % --------------------------------------
        % %% Settings
        % --------------------------------------------
        
        % each design var will control the density and volume fraction
        % material of several clustered elements.
        doUseMultiElePerDV; % 1= true, 0 = false do use multiple elements per (1) design varriable.
        averageMultiElementStrain; % 1= true, 0 = false instead of making 1 large strain field, make sevral and average the sensitivies
        singleMesoDesign; % 1 = yes, 0 = true.
        numXElmPerDV= 2; % Number of elements in the X direction for 1 (per) design varriable.
        numYElmPerDV =2; % Number of elements in the X direction for 1 (per) design varriable.
        numVarsX;
        numVarsY;
        
        mesoplotfrequency=100; % how often to plot the meso level design.
        
        nelxMeso = 5;
        nelyMeso = 5;
        
        nelx = 40; % 40 # of elements in the x direcction, must be a multiple of numXElmPerDV
        nely = 20; % 18 number of elements in the y direction, must be a multiple of numYElmPerDV
        penal = 3; % penality used for the SIMP method
        rmin = 2; % smoothing radius for sensitivity smoothing.
        % Optimization mode and configurations
        mode =4; % 1 = topology only, 2 = material optimization only. 3 = both,4 = meso only. 5 = meso-structure testing
        referenceTemperature = 0; % for thermal expansion, assume that there is not strain when at this temperature.
        addThermalExpansion = 0; % Set to 1 to incorporate thermal expansion
        timestep = 0.1; % time step for the volume fraction update algorithm
        volFractionDamping = 0.1;
        iterationsPerPlot = 5;
        w1 = 1; % weight elastic for multi-objective, % do not set to zero, instead set to 0.0001. 
        w2=0;
        voidMaterialDensityCutOff = 0.3; % everything below this density is considered void.
        
        noNewMesoDesignDensityCutOff = 0.15; % any densities below this will not be redesigned. Having a different value than voidMaterialDensityCutOft helps stabalize the algorithm on the structure edges.
        
        % Plotting information
        doPlotVolFractionDesignVar = 0;
        doPlotTopologyDesignVar = 0;
        doPlotHeat = 0;
        doPlotHeatSensitivityTopology = 0;
        doPlotStress = 0;
        doPlotFinal = 1;
        doPlotMetrics = 0;
        doSaveDesignVarsToCSVFile = 0; % set to 1 to write plotFinal csv file instead
        doPlotAppliedStrain = 0; % For debugging only
        v1 = 0.1; % amount of material 1 to use. default to 10%
        v2 = 0.3; % amount of material 2 to use. default to 30%, reduced so there is less meso structures to compute
        totalVolume; % = v1+v2;
        
        % ----------------
        % Computational settings
        % ------------------
        iterationNum=0; %  used for parallel computing.
        maxFEACalls = 50;
        maxMasterLoops = 30;
        terminationAverageCount = 5; % the average change over this number of iterations must be below the termination criteria
        terminationCriteria = 0.001; % if the normalized average change over  terminationAverageCount of iterations is below this value then termainted. ie. 1% change
        % not much faster.
        useGPU = 0; % set to 1 to try to solve matrix using gpu
        parallel =0; % set to 1 to use parfor while preforming the meso design
        numWorkerProcess = 8;
        useCommandLineArgs = 0;
        
        % -----------------
        % Use different mixture rules for effective elastic properteis
        % 1. Simple linear interpolation, Vigot rule of miztures E = w(E1)*(1-w)*E2
        % 2. Hashin–Shtrikam law (average of upper and lower boundary)
        % 3. Reuss -rule, 1/E = w/E1+(1-w)/E2 (not implemented yet)
        % 4. Mori and Tanaka, metal ceramic composite
        % ---------------------
        elasticMaterialInterpMethod =2;
        % -----------------
        % Use different mixture rules for effective Heat properteis
        % 1. Simple linear interpolation, Vigot rule of miztures E = w(E1)*(1-w)*E2
        % 3. Reuss -rule, 1/E = w/E1+(1-w)/E2 (not implemented yet)
        % 4. Kingery's, metal ceramic composite
        % 5. Hashin–Shtrikam law (average of upper and lower boundary)
        % ---------------------
        heatMaterialInterpMethod = 5;
        
        %          loadingCase = [113]; % left clamped
        %           loadingCase = [111 112 113]; % left clamped
        %            loadingCase = [111 112 ]; % left clamped
        %          loadingCase = [111 120 121]; % up, down, right in top right corrner, left clamp.
        %         loadingCase = [111 120]; % up, down, right in top right corrner, left clamp.
        %              loadingCase = [1];
        
      %  loadingCase = [300 301 302 303 304 305]; % shoe
         loadingCase = [400 401 402 403 404 405]; % bridge
%            loadingCase = [404]; % bridge
        %    loadingCase = [302 305]; % shoe
        
        % --------------
        % Meso tiling info
        %--------------
        numTilesX = 3;
        numTilesY = 3;
        %         sensitivityTile = 5; % use this tile to calcualte the sensitivity
        plotSensitivityWhilerunning = 0;
        
        macro_meso_iteration = 0; % master loop of the whole macro meso system
        
        mesoAddAdjcentCellBoundaries = 0; % global property
        useAjacentLocal = 0; % local property that maybe be turned on or off
        
        
    end
    
    methods
        function obj = Configuration()
            % Constructor.
        end
        
        function [obj]= UpdateVolTargetsAndObjectiveWeights(obj)
            % Do not allow w1 to be zero. Divide by zero messes stuff up.
            if(obj.w1 ==0)
                obj.w1=0.00001;
            end
            obj.w2  = 1- obj.w1; % weight heat transfer
            obj.totalVolume = obj.v1+obj.v2;
        end
        
        function obj = CalculateDesignVarsPerFEAelement(obj)
            
            if(obj.doUseMultiElePerDV ==1)
                if(mod(obj.nely,obj.numYElmPerDV)~=0 || mod(obj.nelx,obj.numXElmPerDV)~=0 )
                    disp('nely and numYElmPerDV or nelx and numXElmPerDV not compatible Exiting MATLAB')
                    exit
                end
                
                obj.numVarsX = obj.nelx/obj.numXElmPerDV;
                obj.numVarsY = obj.nely/obj.numYElmPerDV;
                
            end
        end
        
        % -----------------------------
        % Given design var position, get list of elements X,Y that are
        % controlled by the design var. Since there are multiple elements
        % the x and y returned values are an array.
        %
        % Only applicable when multiple elements per design var is true.
        %
        % -----------------------------
        function [eleXnums, eleYnums,xNodeNums,yNodeNums, macroXdesignVarindex,macroYdesignVarindex] = GetElementsForDesignVar(obj,designvarNumber)
            % 1. Get the x,y position of the design var.
            % 2. multiply by elements per design var to get positions.
            
            numVarsinRow =obj.numVarsX;
            numVarsinColumn =obj.numVarsY;
            ydesignVar =   floor( designvarNumber/numVarsinRow)+1;
            xdesignVar = mod(designvarNumber/numVarsinColumn);
            
            macroXdesignVarindex = xdesignVar;
            macroYdesignVarindex = ydesignVar;
            
            xStart = xdesignVar*obj.numXElmPerDV-(obj.numXElmPerDV);
            yStart = ydesignVar*obj.numYElmPerDV-(obj.numYElmPerDV);
            
            numElementsPerDV = obj.numXElmPerDV*obj.numYElmPerDV;
            eleXnums = zeros(numElementsPerDV,1);
            eleYnums = zeros(numElementsPerDV,1);
            
            count = 1;
            for j = 1:obj.numYElmPerDV
                for i = 1:obj.numXElmPerDV
                    eleXnums(count) = xStart+i;
                    eleYnums(count) = yStart+j;
                    
                    count =count+1;
                end
            end
            
            
            
        end
    end
end