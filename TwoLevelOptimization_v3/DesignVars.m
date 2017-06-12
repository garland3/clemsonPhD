classdef DesignVars
    % Design varriables and support temp varriables are stored in this class.
    
    properties
        % --------------------------------
        % Design Var Arrays
        % --------------------------------
        x; % the "density" at each element
        w; % the volume fraction at each element
        Exx; % the Exx of the orthogonal material
        Eyy; % the Exx of the orthogonal material
        t; %theta, the rotation of the orthogonal material
        E12;
        E33;
        
        % Sub System Copies of design Vars and Lagrangian Multipliers
        lambdaExx=0; % lambda for the Exx
        lambdaEyy=0; % lambda for the Eyy
        lambdaTheta=0; % lambda for the Theta
        penaltyExx=0; % penality value used for the Exx Augmented Lagrangian Multiplier
        penaltyEyy=0; % penality value used for the  Eyy Augmented Lagrangian Multiplier
        penaltyTheta=0; %penality value used for the Theta Augmented Lagrangian Multipliery=
        
        penaltyMesoDensity=0;
        lagrangianMesoDensity=0.1;
        
        ExxSub=0; % Sub system copies of design var
        EyySub=0; % Sub system copies of design var
        thetaSub=0; % Sub system copies of design var
        maxElemStraniEnergy=0;
        
        % ResponseSurfaceCoefficents=[];
        % These are the starting values.
        %  ResponseSurfaceCoefficents=[    -0.0449   1.045e-05   1.045e-05   8.433e-26  -1.045e-10   1.789e-25];
        ResponseSurfaceCoefficents=0;
        averageMesoDensity=0;
        ExxSysAndSubDiffSummed=0;
        EyySysAndSubDiffSummed=0;
        ThetaSysAndSubDiffSummed=0;
        
        % Density off Set (Actual-Predicted)
        densityOffsetArray=0;
        
        
        % Optimization vars (volume fraction optimization)
        lambda1 = 250;
        mu1 = 1;
        c = 0; % objective.
        cCompliance = 0;
        cHeat = 0;
        storeOptimizationVar= []
        
        % --------------------------
        % Support vars
        % -------------------------------
        %         xold; %
        sensitivityElastic; % Sensitivity 1, elastic
        sensitivityElasticPart2; % Needed for E_yy optimization
        sensitivityElasticE12;
        sensitivityElasticE33;
        
        sensitivityHeat; % Sensitivity 2, heat
        currentVol1Fraction;
        currentVol2Fraction;
        
        actualAverageE=0;
        
        %complianceSensitivity; %
        totalStress;
        dc; % Derivative of c (hence dc). C is the objective.
        dcSum; % derivitive of meso element sensitivies when designing a single meso structure for the whole macro structure.
        %g1elastic; % Derivative of c with respect to a material change for the elastic
        %g1heat; %  Derivative of cHeat with respect to a material change for the heat
        IEN; % element to node map. Save this matrix, so it does not need to be recalculated every time.
        XLocations; %=zeros(numNodesInRow,numNodesInColumn);
        YLocations; %=zeros(numNodesInRow,numNodesInColumn);
        globalPosition; %  = zeros(nn,2);
        NodeToXYArrayMap; % map of node numbers to their X,Y position in FEA arrays
        elementXYposition; % or othe position in the X ad W matrixes
        
        U_heatColumn; % temperature matrix, gives the temperature at each node (not at the element centers)
        U; % displacement matrix
        B; %  holding the displacement strain matrix
        
        maxF = 0;
        maxU = 0;
        
        
        % meso informatoin
        nelxTile;
        nelyTile
        xTile;
        wTile;
        IENTile;
        XLocationsTile;
        YLocationsTile;
        globalPositionTile;
        NodeToXYArrayMapTile;
        UTile;
        d11;
        d12;
        d22;
        d33;
        De11
        De12
        De22
        De33;
        
        mesoStructNTCmask;
        mesoFEACalls=0;
        
        
        % When multiple elements are controlled by a single design var,
        % this array gives you the design var number.
        % Give it a X,Y position and it will tell you the the rows are the
        % element numbers.
        elementToDesignVarMap;
    end
    
    methods
        % Constructur method
        function obj = DesignVars(config)
            
            % Get the B matrix, the E,v,G do not matter and are not used
            % in the B calculation, so set them to 1.
            E = 1; v= 1; G = 1; strain = [];
            [~, ~, B_out] = elK_elastic(E,v, G,strain,[]);
            obj.B = B_out;
            
            macroDesignMode = 0;
            if(config.mode<100 || config.mode==111 || config.mode==112) % mode 111 and  112 is meso validation target generation
                macroDesignMode=1;
            end
            
            if ( macroDesignMode==1)
                % ------------------------------
                % Initialize Arrays used for Macro Optimization
                % ------------------------------
                obj.x(1:config.nely,1:config.nelx) =config.v1+config.v2; % artificial density of the elements
                obj.w(1:config.nely,1:config.nelx)  =config.v1; % actual volume fraction composition of each element
                %                   obj.w(1:config.nely,1:config.nelx)  =ones(config.nely,config.nelx);
                %                 obj.d(1:config.nely,1:config.nelx)  =ones(config.nely,config.nelx)*0.8; % orthotropic masterial distribution
                obj.Exx(1:config.nely,1:config.nelx)  =ones(config.nely,config.nelx);
                obj.Eyy(1:config.nely,1:config.nelx)  =ones(config.nely,config.nelx);
                obj.t(1:config.nely,1:config.nelx)  =ones(config.nely,config.nelx)*0; % rotation of the orthotropic material
                
                obj.sensitivityElastic(1:config.nely,1:config.nelx) = 0;
                obj.sensitivityElasticPart2 (1:config.nely,1:config.nelx) = 0;
                obj.sensitivityHeat(1:config.nely,1:config.nelx) = 0;
                if (config.doPlotStress == 1)
                    obj.totalStress(1:config.nely,1:config.nelx) = 0;
                end
                obj=obj.CalcIENmatrix(config);
                obj=obj.CalcNodeLocation(config);
                obj=obj.PreCalculateXYmapToNodeNumber(config);
                obj=obj.CalcElementXYposition(config);
                %     matProp=  matProp.ReadConstitutiveMatrixesFromFiles(config);
                
                
                % Start them out as 0. Set later if iteration 2 or higher
                obj.lambdaExx= obj.Exx*0; % lambda for the Exx
                obj.lambdaEyy= obj.Exx*0; % lambda for the Eyy
                obj.lambdaTheta=obj.Exx*0; % lambda for the Theta
                obj.penaltyExx= obj.Exx*0; % penality value used for the Exx Augmented Lagrangian Multiplier
                obj.penaltyEyy= obj.Exx*0; % penality value used for the  Eyy Augmented Lagrangian Multiplier
                obj.penaltyTheta= obj.Exx*0; %penality value used for the Theta Augmented Lagrangian Multiplier
                
                obj.ExxSub=obj.Exx*0; % Sub system copies of design var
                obj.EyySub=obj.Exx*0; % Sub system copies of design var
                obj.thetaSub=obj.Exx*0; % Sub system copies of design var
                obj.densityOffsetArray=obj.Exx*0;
                
                if(config.useTargetMesoDensity==1)
                    obj.Exx = ones(size(obj.Exx))*100000*config.targetExxEyyDensity;
                    %      DV.Exx = ones(size(DV.Exx))*3.277340e+04;
                else
                    obj.Exx = ones(size(obj.Exx))*config.targetAvgExxEyy;
                end
                obj.Eyy = obj.Exx ;
                
                
            end
            
            if(config.anisotropicMat==1)
                obj.E12 = obj.Exx;
                obj.E33 = obj.Exx;
                obj. sensitivityElasticE12= obj.sensitivityElastic;
                obj.sensitivityElasticE33= obj.sensitivityElastic;
            end
            
            if(config.useThetaInSurfaceFit==1)
                %                 obj. ResponseSurfaceCoefficents=[  9.99999999779786e-06 9.9999086180473e-06 9.99990080613124e-06 1.00000000554909e-05 -1.33226777965037e-10 -1.94015928184367e-10 9.99999972290092e-06 4.71726205107711e-11 9.99990147075606e-06 9.99982147731046e-06];
                
                obj. ResponseSurfaceCoefficents=[ -0.1357    0.3071    1.1000   -0.1042    0.5286   -0.4000    0.1327   -0.4000   -0.1273    0.1273];
            else
                % obj. ResponseSurfaceCoefficents=[ 1.0000000000463e-05 9.99988184437107e-06 9.9998491550433e-06 -3.40115537230351e-11 -5.52110060132392e-12 -3.81038581303971e-11];
                %  obj.ResponseSurfaceCoefficents=[    -0.0449   1.045e-05   1.045e-05   8.433e-26  -1.045e-10   1.789e-25];
                
                %Data from the estimated values that I came up with 
%                 obj.ResponseSurfaceCoefficents=[     -0.0449    1.0449    1.0449    0.0000   -1.0449    0.0000];

            % from lookup table using the validation data
            % fit is a funciton of Exx/matProp.E_material1 and Eyy/matProp.E_material1
            % poly33 is the function. 
             obj. ResponseSurfaceCoefficents=[    0.1137      1.376    1.353    -0.6908    -1.988      -0.6529      0.1601   0.6243      0.5984   0.1471 ];
             
          %   obj. ResponseSurfaceCoefficents=[     0.0847            1.478           1.452       -0.8023           -2.179           -0.7586            0.2135          0.6952            0.6678            0.1974 ];
            end
            
            
            
        end
        
        % -----------------------------
        % Get the state of the macro optimization saved in csv files.
        %
        % Also, read Exx,Eyy,Theta, Lambdas and ExxSub, EyySub, ThetaSub
        % -----------------------------
        function obj = GetMacroStateVarsFromCSVFiles(obj, config,matProp)
            
            if(config.macro_meso_iteration>1 ||config.mode==90)
                
                if (config.mode==60 ||config.mode==90  )
                    % -------------------------
                    % Read from csv
                    % x (topology var)
                    % Exx System
                    % Eyy System
                    % t (Theta) System
                    %
                    % lambdaExx Array
                    % lambdaEyy Array
                    % lambdaTheta Array
                    %
                    % Exx SubSystem
                    % Eyy SubSystem
                    % t (Theta) SubSystem
                    %
                    % penality Exx
                    % penalty Eyy
                    % penalty theta
                    %
                    % See function SaveMacroProblemStateToCSV and GenerateRhoFunctionOfExxEyy
                    % -------------------------
                    folderNum = config.iterationNum;
                    oldIteration = config.macro_meso_iteration-1;
                    if(config.mode==90)
                         oldIteration = config.macro_meso_iteration;
                    end
                    
                    outnameX = sprintf('./out%i/SIMPdensityfield%i.csv',folderNum,oldIteration);
                    outnameExx = sprintf('./out%i/ExxValues%i.csv',folderNum,oldIteration);
                    outnameEyy = sprintf('./out%i/EyyValues%i.csv',folderNum,oldIteration);
                    outnameTheta = sprintf('./out%i/ThetaValues%i.csv',folderNum,oldIteration);
                    outnamelambdaExx = sprintf('./out%i/lambdaExx%i.csv',folderNum,oldIteration);
                    outnamelambdaEyy = sprintf('./out%i/lambdaEyy%i.csv',folderNum,oldIteration);
                    outnamelambdaTheta = sprintf('./out%i/lambdaTheta%i.csv',folderNum,oldIteration);
                    outnamepenaltyExx = sprintf('./out%i/penaltyExx%i.csv',folderNum,oldIteration);
                    outnamepenaltyEyy = sprintf('./out%i/penaltyEyy%i.csv',folderNum,oldIteration);
                    outnamepenaltyTheta = sprintf('./out%i/penaltyTheta%i.csv',folderNum,oldIteration);
                    outnameExxSubSysValues = sprintf('./out%i/ExxSubSysValues%i.csv',folderNum,oldIteration);
                    outnameEyySubSysValues = sprintf('./out%i/EyySubSysValues%i.csv',folderNum,oldIteration);
%                     outnameThetaSubSysValues = sprintf('./out%i/ThetaSubSysValues%i.csv',folderNum,oldIteration);
                    outnameThetaSubSysValues = sprintf('./out%i/ThetaSubSysValues%i.csv',folderNum,oldIteration);
                     outnameMesoDensities = sprintf('./out%i/densityUsedSubSysValues%i.csv',folderNum,oldIteration);
                  
                    
                    
                    obj.x = csvread(outnameX); % the "density" at each element
                    obj.Exx = csvread(outnameExx); % the Exx of the orthogonal material
                    obj.Eyy = csvread(outnameEyy); % the Eyy of the orthogonal material
                    obj.t = csvread(outnameTheta); %theta, the rotation of the orthogonal material
                    
                    % Sub System Copies of design Vars and Lagrangian Multipliers
                    % starting on iteration 3, get the saved values.
                    % these values are initilized on iteration 2
                    if(config.macro_meso_iteration>2)
                        obj.lambdaExx= csvread(outnamelambdaExx); % lambda for the Exx
                        obj.lambdaEyy= csvread(outnamelambdaEyy); % lambda for the Eyy
                        obj.lambdaTheta= csvread(outnamelambdaTheta); % lambda for the Theta
                        obj.penaltyExx= csvread(outnamepenaltyExx); % penality value used for the Exx Augmented Lagrangian Multiplier
                        obj.penaltyEyy= csvread(outnamepenaltyEyy); % penality value used for the  Eyy Augmented Lagrangian Multiplier
                        obj.penaltyTheta= csvread(outnamepenaltyTheta); %penality value used for the Theta Augmented Lagrangian Multiplier
                    elseif(config.macro_meso_iteration==2)
                        % set the lagrangians to zero for now. They will be
                        % updated later after the penalty is calculated.
                        obj.lambdaExx=zeros(config.nely,config.nelx);
                        obj.lambdaEyy= obj.lambdaExx;
                        obj.lambdaTheta= obj.lambdaExx;
                    end
                    
                    obj.ExxSub=csvread(outnameExxSubSysValues); % Sub system copies of design var
                    obj.EyySub=csvread(outnameEyySubSysValues); % Sub system copies of design var
                    obj.thetaSub=csvread(outnameThetaSubSysValues); % Sub system copies of design var
                    
                    actualDensities=csvread(outnameMesoDensities); % Sub system copies of design var
                    if(config.useTargetMesoDensity==1&& config.mode~=90)
                        o = Optimizer;
                        p = plotResults;
                        figure
                        [~, ~,predictedDensities] = o.CalculateDensitySensitivityandRho(obj.ExxSub/matProp.E_material1,obj.EyySub/matProp.E_material1,obj.thetaSub,obj.x ,obj.ResponseSurfaceCoefficents,config,matProp,obj.densityOffsetArray);
                        subplot(3,1,1);
                        p.PlotArrayGeneric(actualDensities,'Actual Meso Densities');
                        subplot(3,1,2);
                        p.PlotArrayGeneric(predictedDensities,'Predicted Meso Densities');
                        subplot(3,1,3);
                        p.PlotArrayGeneric(actualDensities-predictedDensities,'Actual - Predicted Meso Densities');
                        nameGraph = sprintf('./PredictedVSActualMesoDensities%i.png', config.macro_meso_iteration);
                        print(nameGraph,'-dpng');
                        
                        obj.densityOffsetArray=(actualDensities-predictedDensities)*0.5;
                        close all
                        p.PlotArrayGeneric( obj.densityOffsetArray,'Offset Densities');
                        nameGraph = sprintf('./Predicted_OffsetDensitiesForIteration%i.png', config.macro_meso_iteration);
                        print(nameGraph,'-dpng');
                        %                             obj.densityOffsetArray=actualDensities*0;
                    end
                    
                    %
                    %                     if(config.macro_meso_iteration>1)
                    %                         nameArray = sprintf('./out%i/ExxEyyRhoFitCoefficients%i.csv',folderNum, oldIteration);
                    %                          xxx= dlmread(nameArray);
                    %                          sprintf('%.12f',xxx(1));
                    %                         obj.ResponseSurfaceCoefficents=xxx;
                    %                     end
                    
                    
                elseif(config.mode==50)
                    % -------------------------
                    % (old) if Volume Fraction optimization
                    % -------------------------
                    folderNum = config.iterationNum;
                    previousIterationNum = config.macro_meso_iteration-1;
                    
                    % get the topology densities.
                    outname = sprintf('./out%i/densityfield%i.csv',folderNum,previousIterationNum);
                    obj.x = csvread(outname);
                    
                    % get the volume fraction optimization vars
                    outname = sprintf('./out%i/volfractionfield%i.csv',folderNum,previousIterationNum);
                    obj.w=csvread(outname);
                    
                    % get the lambda value
                    outname = sprintf('./out%i/lambda%i.csv',folderNum,previousIterationNum);
                    obj.lambda1=csvread(outname);
                elseif(config.mode==1 && config.multiscaleMethodCompare==1)
                     folderNum = config.iterationNum;
                    oldIteration = config.macro_meso_iteration-1;                    
                    outnameX = sprintf('./out%i/SIMPdensityfield%i.csv',folderNum,oldIteration);
                    obj.x = csvread(outnameX); % the "density" at each element
                    
                    
                else
                   fprintf('First iteration, no macro state vars to read');
                end
            end
        end
        
        
        
        % -----------------------------
        % UpdatePenaltyAndLagrangianValues
        % -----------------------------
        function obj = UpdatePenaltyAndLagrangianValues(obj, config,matProp)
            
            if(config.mode==90)
                return
            end
            
            ne = config.nelx*config.nely;
            diffExx = obj.ExxSub-obj.Exx;
            diffEyy = obj.EyySub-obj.Eyy  ;
            diffTheta = obj.thetaSub-obj.t;
            %              diffExx = -obj.Exx+obj.ExxSub;
            %             diffEyy = -obj.Eyy + obj.EyySub;
            %             diffTheta = -obj.t+obj.thetaSub;
            
            
            
            % if the first time, calculate the initial penalty values.
            if(config.macro_meso_iteration>=2)
                % Get displacement fields
                oldIteration = config.macro_meso_iteration-1;
                folderNum = config.iterationNum;
                outname = sprintf('./out%i/displacement%i.csv',folderNum,oldIteration);
                UpreviousIteration =  csvread(outname);
                
                % read the Exx sensitivity field
                outname = sprintf('./out%i/sensitivityElastic%i.csv',folderNum,oldIteration);
                ExxSensitivity= csvread(outname);
                
                % read the Eyy sensitivity field
                outname = sprintf('./out%i/sensitivityElasticPart2%i.csv',folderNum,oldIteration);
                EyySensitivity = csvread(outname);
                
%             end
%             if(config.macro_meso_iteration>=2)
                
                % loop over the elements.
                
                for e = 1:ne
                    
                    [xPos,yPos]= obj.GivenNodeNumberGetXY(e);
                    
                    nodes1 = obj.IEN(e,:);
                    xNodes = nodes1*2-1;
                    yNodes = nodes1*2;
                    
                    % I cannot use the union, or else the order get messed up. The order
                    % is important. Same in the actual topology code when you are
                    % calculating the objectiv
                    NodeNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
                    
                    E12_local = 1;
                    E33_local = 1;
                    KE = matProp.getKMatrixTopExxYyyRotVars(config,obj.x(yPos,xPos),obj.Exx(yPos,xPos), obj.Eyy(yPos,xPos),obj.t(yPos,xPos),obj.w(yPos,xPos),E12_local,E33_local,e);
                    
                    
                    % get the displacements and calculate the strain energy
                    strainEnergy = 0;
                    [~, t2] = size(config.loadingCase);
                    for loadcaseIndex = 1:t2
                        u_local = transpose(UpreviousIteration(loadcaseIndex,NodeNumbers));
                        strainEnergy_temp=transpose(u_local)*KE*u_local;
                        strainEnergy=strainEnergy+strainEnergy_temp;
                        
                        
                    end
                    
                    % ---------------------------------------------
                    % Calculate the penalty values!!
                    % ---------------------------------------------
                    omegaLocal =config.Omega;
                    omegaLocal=omegaLocal*0.5;
                    
                    if(config.macro_meso_iteration>=3)
                        multiplier = 4*(config.macro_meso_iteration-2);
                        omegaLocal=omegaLocal*multiplier;
                             omegaLocal = min(omegaLocal,2);
                    end
                    %                     obj.penaltyExx(yPos,xPos)=min(2*omegaLocal*strainEnergy/ abs( diffExx(yPos,xPos)),strainEnergy);
                    %                     obj.penaltyEyy(yPos,xPos)=min(2*omegaLocal*strainEnergy/ abs (diffEyy(yPos,xPos)),strainEnergy);
                    %                     obj.penaltyTheta(yPos,xPos)=min(2*omegaLocal*strainEnergy/ abs( diffTheta(yPos,xPos)),strainEnergy); % max penalty of strainEnergy*100
                    %
                    %                     obj.penaltyExx(yPos,xPos)=omegaLocal*ExxSensitivity(yPos,xPos);
                    %                     obj.penaltyEyy(yPos,xPos)=omegaLocal*EyySensitivity(yPos,xPos);
                    obj.penaltyTheta(yPos,xPos)=omegaLocal*strainEnergy;
                    %
                end
                
            end
            
            %             maxP = max(max(max( obj.penaltyExx)),max(max(obj.penaltyEyy)));
            %               obj.penaltyExx=config.Omega*  obj.penaltyExx/maxP;% divide by largest, scale by the Omega
            %                 obj.penaltyEyy= config.Omega* obj.penaltyEyy/maxP;
            %
            %                 maxPTheta = max(max( obj.penaltyTheta));
            %                  obj.penaltyTheta=config.Omega* obj.penaltyTheta/maxPTheta;
            
            
            
            % ---------------------------------------------
            % Calculate the penalty values!! if on iteration 3 or
            % greater
            % ---------------------------------------------
%             if(config.macro_meso_iteration>2)
%                 
%                 updateMultiplier = 4;
%                 for e = 1:ne
%                     [xPos,yPos]= obj.GivenNodeNumberGetXY(e);
%                     %                     obj.penaltyExx(yPos,xPos)=   obj.penaltyExx(yPos,xPos)*updateMultiplier;
%                     %                     obj.penaltyEyy(yPos,xPos)= obj.penaltyEyy(yPos,xPos)*updateMultiplier;
%                     obj.penaltyTheta(yPos,xPos)=  obj.penaltyTheta(yPos,xPos)*updateMultiplier;
%                 end
%                 
%             end
            
            
            % ---------------------------------------------
            % Calculate the lagrangian values.
            % ---------------------------------------------
            obj.lambdaExx= zeros(size(obj.Exx));
            obj.lambdaEyy=  zeros(size(obj.Exx));
            
            %             for e = 1:ne
            %                 [xPos,yPos]= obj.GivenNodeNumberGetXY(e);
            %                 deltaT = 0.1;
            % %                 obj.lambdaExx(yPos,xPos)= obj.lambdaExx(yPos,xPos)+deltaT *diffExx(yPos,xPos);
            % %                 obj.lambdaEyy(yPos,xPos)=  obj.lambdaEyy(yPos,xPos)+deltaT*diffEyy(yPos,xPos);
            %
            %                 obj.lambdaTheta(yPos,xPos)=   obj.lambdaTheta(yPos,xPos)+obj.penaltyTheta(yPos,xPos) *diffTheta(yPos,xPos);
            %
            %
            %
            % %                 obj.lambdaExx(yPos,xPos)=0; %obj.lambdaExx(yPos,xPos)+obj.penaltyExx(yPos,xPos) *diffExx(yPos,xPos);
            % %                 obj.lambdaEyy(yPos,xPos)=  0;%obj.lambdaEyy(yPos,xPos)+obj.penaltyEyy(yPos,xPos) *diffEyy(yPos,xPos);
            % %                 obj.lambdaTheta(yPos,xPos)=  0;% obj.lambdaTheta(yPos,xPos)+obj.penaltyTheta(yPos,xPos) *diffTheta(yPos,xPos)
            %             end
            
            
            % calculate the penalties now
            if(config.macro_meso_iteration>=2)
                deltaT=1;
                obj.lambdaExx=diffExx*deltaT;
                obj.lambdaEyy=diffEyy*deltaT;
                
                smallestLambdExx = min(min(obj.lambdaExx));
                smallestLambdEyy = min(min(obj.lambdaEyy));
                smallestOfTwo = min(smallestLambdExx,smallestLambdEyy)-1;
                omegaLocal =config.Omega;
                if(config.macro_meso_iteration>=3)
                    multiplier = 2*(config.macro_meso_iteration-2);
                 
                    omegaLocal=omegaLocal*multiplier;
                       omegaLocal = min(omegaLocal,2);
                end
                
                obj.penaltyExx=abs(omegaLocal*ExxSensitivity./(obj.lambdaExx-smallestOfTwo));
                obj.penaltyEyy=abs(omegaLocal*EyySensitivity./(obj.lambdaEyy-smallestOfTwo));
                
                %                 obj.penaltyExx=obj.penaltyExx*25;
                %                 obj.penaltyEyy=obj.penaltyEyy*25;
                
                % d =matProp.E_material1
                %                 p = ones(size(EyySensitivity))*max(max(ExxSensitivity))/matProp.E_material1;
                % %                  obj.penaltyExx=abs(omegaLocal*ExxSensitivity./(obj.lambdaExx-smallestOfTwo));
                % %                 obj.penaltyEyy=abs(omegaLocal*EyySensitivity./(obj.lambdaEyy-smallestOfTwo));
                %                 obj.penaltyExx=p;
                %                 obj.penaltyEyy=p;
                
            end
            
            if 1==0
                figure(1)
                p = plotResults;
                figure
                subplot(2,3,1)
                p.PlotArrayGeneric( obj.penaltyExx, ' obj.penaltyExx')
                subplot(2,3,2)
                p.PlotArrayGeneric( obj.penaltyEyy, ' obj.penaltyEyy')
                subplot(2,3,3)
                p.PlotArrayGeneric( obj.penaltyTheta, ' obj.penaltyTheta')
                
                subplot(2,3,4)
                p.PlotArrayGeneric( obj.lambdaExx, ' obj.lambdaExx')
                subplot(2,3,5)
                p.PlotArrayGeneric( obj.lambdaEyy, ' obj.lambdaEyy')
                subplot(2,3,6)
                p.PlotArrayGeneric( obj.lambdaTheta, ' obj.lambdaTheta')
                drawnow
            end
        end
        
        
        
        % IEN holds the node numbers for each element.
        % Each row is a new element
        % The first column is element 1's global node number
        % Second column is elemnt 2's global node number
        %  and ....
        %
        % Element nodes are as follows
        %
        %  4 ---- 3
        %  |      |
        %  |      |
        %  1 ---- 2
        %
        %
        % Let the x direction be the first dof and y direction the second
        % dof
        %
        %   y = 2
        %   /\
        %   |
        %   |
        %   *----> x = 1
        %
        %
        %
        % ---------------------------------------------
        % Global matrix nodes are as follows
        % row = nelx+1
        % col = nely+1
        %
        % col*row+1-col*row+2-col*row+3-col*row+4...col*row+row-1-col*row+row
        % .            .         .         .             .            .
        % .            .         .         .             .            .
        % .            .         .         .             .            .
        % |            |         |         |             |            |
        % 2*row+1 - 2*row+2 - 2*row+3 - 2*row+4 ... 2*row+row-1 - 2*row+row
        % |            |         |         |             |            |
        % 1*row+1 - 1*row+2 - 1*row+3 - 1*row+4 ... 1*row+row-1 - 1*row+row
        % |            |         |         |             |            |
        % 0*row+1 - 0*row+2 - 0*row+3 - 0*row+4 ... 0*row+row-1 - 0*row+row
        function obj =  CalcIENmatrix(obj,config)
            
            count = 1;
            elementsInRow = config.nelx+1; % think this is actually "nodes in a row. "
            nn = (config.nelx+1)*(config.nely+1); % number of nodes
            obj.IEN = zeros(nn,4);
            % Each row, so nely # of row
            for i = 1:config.nely
                rowMultiplier = i-1;
                % Each column, so nelx # of row
                for j= 1:config.nelx
                    obj.IEN(count,:)=[rowMultiplier*elementsInRow+j, ...
                        rowMultiplier*elementsInRow+j+1, ...
                        (rowMultiplier +1)*elementsInRow+j+1,...
                        (rowMultiplier +1)*elementsInRow+j];
                    count = count+1;
                end
            end
        end
        
        
        % Calcualte the node locations and put the information
        % into an array. Needed for the FEA
        % Calculate it here, so it only need to be calculated once.
        %
        % Think these are actually node locations
        function obj = CalcNodeLocation(obj,config)
            nn = (config.nelx+1)*(config.nely+1); % number of nodes
            
            numNodesInRow = config.nelx + 1;
            numNodesInColumn = config.nely + 1;
            %             obj.XLocations=zeros(numNodesInRow,numNodesInColumn);
            %             obj.YLocations=zeros(numNodesInRow,numNodesInColumn);
            
            obj.globalPosition = zeros(nn,2);
            count = 1;
            for i = 1:numNodesInColumn  % y
                for j= 1:numNodesInRow % x
                    obj.globalPosition(count,:) = [j-1 i-1];
                    count = count+1;
                    %                     obj.XLocations(j,i) = j-1;
                    %                     obj.YLocations(j,i) = i-1;
                end
            end
        end
        
        % Pre Calculate a map of the XY array coordinates for each node
        % number.
        function obj= PreCalculateXYmapToNodeNumber(obj ,config)
            nn = (config.nelx+1)*(config.nely+1); % number of nodes
            obj.NodeToXYArrayMap = zeros(nn,2);
            count = 1;
            for i = 1:config.nely
                for j= 1:config.nelx
                    obj.NodeToXYArrayMap(count,:) = [i,j];
                    count = count+1;
                end
            end
        end
        
        % -----------------------------
        %
        % Get the design var position when given an element X,Y position.
        %
        % Only applicable when multiple elements per design var is true.
        %
        % -----------------------------
        %         function [x,y] = GetDesignVarPositionGivenXYElement(obj,config,xelm,yelm)
        %             x = ceil(xelm/config.numXElmPerDV);
        %             y = ceil(yelm/config.numYElmPerDV);
        %         end
        
        
        
        % -----------------------------
        % Calcualte the Center of each element and put the information
        % into an array. Needed for the FEA
        % Calculate it here, so it only need to be calculated once.
        % -----------------------------
        function obj = CalcNodeLocationMeso(obj,config)
            nn = (config.nelx)*(config.nely); % number of nodes
            
            numNodesInRow = config.nelx ;
            numNodesInColumn = config.nely ;
            obj.XLocations=zeros(numNodesInRow,numNodesInColumn);
            obj.YLocations=zeros(numNodesInRow,numNodesInColumn);
            
            obj.globalPosition = zeros(nn,2);
            count = 1;
            for i = 1:(numNodesInColumn)  % y
                for j= 1:(numNodesInRow) % x
                    obj.globalPosition(count,:) = [j-1 i-1];
                    count = count+1;
                    obj.XLocations(j,i) = j-1;
                    obj.YLocations(j,i) = i-1;
                end
            end
        end
        
        % -----------------------------
        % -----------------------------
        function obj = CalcNodeLocationMeso_Tile(obj,config)
            nn = (obj.nelxTile)*(obj.nelyTile); % number of nodes same as number of elements since it wraps
            %             nn = (config.nelx)*(config.nely); % number of nodes
            
            numNodesInRow = obj.nelxTile ;
            numNodesInColumn = obj.nelyTile;
            obj.XLocations=zeros(numNodesInRow,numNodesInColumn);
            obj.YLocations=zeros(numNodesInRow,numNodesInColumn);
            
            obj.globalPosition = zeros(nn,2);
            count = 1;
            for i = 1:(numNodesInColumn)  % y
                for j= 1:(numNodesInRow) % x
                    obj.globalPositionTile(count,:) = [j-1 i-1];
                    count = count+1;
                    obj.XLocationsTile(j,i) = j-1;
                    obj.YLocationsTile(j,i) = i-1;
                end
            end
        end
        
        
        
        % ---------------------------
        %
        %  Calcualte the element node map (IEN) for the homogenization meso
        %  structure case. BAsically the map needs to have the displacement
        %  field be periodic, so ti must loop back on itself on the edges.
        %
        % ---------------------------
        function obj =  CalcElementNodeMapmatrixWithPeriodicXandY(obj,config)
            
            count = 1;
            elementsInRow = config.nelx;
            nn = (config.nelx)*(config.nely); % number of nodes same as number of elements since it wraps
            obj.IEN = zeros(nn,4);
            % Each row, so nely # of row
            for i = 1:config.nely
                rowMultiplier = i-1;
                % Each column, so nelx # of row
                for j= 1:config.nelx
                    
                    % normal case
                    if(j ~= config.nelx && i ~= config.nely )
                        obj.IEN(count,:)=[rowMultiplier*elementsInRow+j, ...
                            rowMultiplier*elementsInRow+j+1, ...
                            (rowMultiplier +1)*elementsInRow+j+1,...
                            (rowMultiplier +1)*elementsInRow+j];
                    elseif(j == config.nelx && i ~= config.nely )
                        % On the right side of the mesh case
                        
                        % loop back arround case case
                        obj.IEN(count,:)=[rowMultiplier*elementsInRow+j, ...
                            rowMultiplier*elementsInRow+1, ... % note the difference
                            (rowMultiplier +1)*elementsInRow+1,... % note the difference
                            (rowMultiplier +1)*elementsInRow+j];
                    elseif(j ~= config.nelx && i == config.nely )
                        % On the top side of the mesh case
                        
                        obj.IEN(count,:)=[rowMultiplier*elementsInRow+j, ...
                            rowMultiplier*elementsInRow+j+1, ...
                            0+j+1,... % note the difference
                            0+j]; % note the difference
                        
                    elseif(j == config.nelx && i == config.nely )
                        % On the top right  side of the mesh case
                        obj.IEN(count,:)=[rowMultiplier*elementsInRow+j, ...
                            rowMultiplier*elementsInRow+1, ...
                            1,... % back to the node 1
                            j];
                        
                    end
                    count = count+1;
                end
            end
        end
        
        
        
        
        % --------------------------------------------
        %
        %       ELEMENT TO NODE MAP FOR TILED PERIODIC MESO
        %
        % --------------------------------------------
        function obj =  CalcElementNodeMapmatrixWithPeriodicXandY_Tile(obj,config)
            
            count = 1;
            elementsInRow=obj.nelxTile;% elementsInRow = config.nelx;
            nn = (obj.nelxTile)*(obj.nelyTile); % number of nodes same as number of elements since it wraps
            obj.IENTile = zeros(nn,4);
            % Each row, so nely # of row
            for i = 1:obj.nelyTile
                rowMultiplier = i-1;
                % Each column, so nelx # of row
                for j= 1:obj.nelxTile
                    
                    % normal case
                    if(j ~= obj.nelxTile && i ~= obj.nelyTile )
                        obj.IENTile(count,:)=[rowMultiplier*elementsInRow+j, ...
                            rowMultiplier*elementsInRow+j+1, ...
                            (rowMultiplier +1)*elementsInRow+j+1,...
                            (rowMultiplier +1)*elementsInRow+j];
                    elseif(j == obj.nelxTile && i ~= obj.nelyTile )
                        % On the right side of the mesh case
                        
                        % loop back arround case case
                        obj.IENTile(count,:)=[rowMultiplier*elementsInRow+j, ...
                            rowMultiplier*elementsInRow+1, ... % note the difference
                            (rowMultiplier +1)*elementsInRow+1,... % note the difference
                            (rowMultiplier +1)*elementsInRow+j];
                    elseif(j ~= obj.nelxTile && i == obj.nelyTile )
                        % On the top side of the mesh case
                        
                        obj.IENTile(count,:)=[rowMultiplier*elementsInRow+j, ...
                            rowMultiplier*elementsInRow+j+1, ...
                            0+j+1,... % note the difference
                            0+j]; % note the difference
                        
                    elseif(j == obj.nelxTile && i == obj.nelyTile )
                        % On the top right  side of the mesh case
                        obj.IENTile(count,:)=[rowMultiplier*elementsInRow+j, ...
                            rowMultiplier*elementsInRow+1, ...
                            1,... % back to the node 1
                            j];
                        
                    end
                    count = count+1;
                end
            end
        end
        
        % Given an X, find the node number
        function number = GetNodeNumberGivenXY(obj,config, x,y)
            numNodesInRow = config.nelx+1;
            % numNodesInColumn = obj.nely+1;
            rowMultiplier = y-1;
            number = rowMultiplier*numNodesInRow+x;
        end
        
        % Given a node number, find the X, Y position (not the physical
        % position, but the matrix location position)
        function [x , y ]= GivenNodeNumberGetXY(obj, nodeNum)
            [result] =  obj.NodeToXYArrayMap(nodeNum,:);
            y = result(1);
            x = result(2);
        end
        
        % Calculates the volume
        function [obj] = CalculateVolumeFractions(obj, config,matProp)
            
            
            %             Xtemp = obj.x;
            %             Xtemp(obj.x>config.voidMaterialDensityCutOff)=1;
            %             Xtemp(obj.x<=config.voidMaterialDensityCutOff)=0;
            %             ne = config.nelx*config.nely;
            %                ne = config.nelx*config.nely;
            %                 averageElasticLocal = (sum(sum(obj.Exx.*obj.x))+sum(sum(obj.Eyy.*obj.x)))/ne;
            %                 volume1 = sum(sum(obj.x))/ne;
            %                 volume2 = averageElasticLocal/matProp.E_material1;
            %                   volume1 = (sum(sum(obj.Exx.*obj.x))+sum(sum(obj.Eyy.*obj.x)))/ne;
            
            ne = config.nelx*config.nely;
            %             neSolid = config.nelx*config.nely*(config.v1+config.v2);
            totalMaterial = sum(sum(obj.x));
%             obj.targetAverageE=(config.v1*matProp.E_material1+config.v2*matProp.E_material2)/(config.v1+config.v2);
            
            if(config.useExxEyy==1)
                %                    avg= 0.5*(obj.Exx+obj.Eyy);
                totalExx =obj.x.*obj.Exx;
                totalEyy = obj.x.* obj.Eyy;
                avgE = (totalExx+totalEyy)/2;
                if(config.testingVerGradMaterail ==1)
                    minE = matProp.E_material2;
                else
                    minE = matProp.E_material2/2;
                end
                obj.w = avgE/(matProp.E_material1-minE);
                %                   for i = 1:config.nelx
                %                       for j = 1:config.nely
                %                            obj.w(j,i) =max(obj.Exx(j,i),obj.Eyy(j,i))/(matProp.E_material1-minE);
                %                       end
                %                   end
                %                   obj.w =max(obj.Exx,obj.Eyy)/matProp.E_material1;
                
                
                obj.actualAverageE= sum(sum(avgE))/totalMaterial;
                
                % volume fraction is over the domain, while target E
                obj.  currentVol1Fraction =sum(sum( obj.x.*obj.w))/ne;
                obj.   currentVol2Fraction =sum(sum( obj.x.*(1-obj.w)))/ne;
                
                
                
            else
                
                
                
                totalMat1 =sum(sum( obj.x.*obj.w*matProp.E_material1));
                totalMat2 =sum(sum( obj.x.*(1-obj.w)*matProp.E_material2));
                % obj.actualAverageE= obj.currentVol1Fraction*matProp.E_material1+  obj. currentVol2Fraction*matProp.E_material2;
                obj.actualAverageE= (totalMat1+totalMat2)/totalMaterial;
                obj.  currentVol1Fraction =sum(sum( obj.x.*obj.w))/ne;
                obj.   currentVol2Fraction =sum(sum( obj.x.*(1-obj.w)))/ne;
            end
            
            
        end
        
        
      
        %%
        
        % --------------------------------------------
        %
        %      Calculate the element XY position.
        %
        % --------------------------------------------
        function obj = CalcElementXYposition(obj,config)
            nelm = config.nelx*config.nely;
            
            obj.elementXYposition =   zeros(nelm,2);
            
            e_count = 1;
            for ely = 1:config.nely
                %                 rowMultiplier = ely-1;
                for elx = 1:config.nelx
                    obj.elementXYposition(e_count,:)=[ely,elx];
                    e_count = e_count+1;
                end
            end
        end
        
        
        % ----------------------------------
        %
        %  RUN THE FEA ANAYLSIS FOR EACH LOADING CASE
        %
        % --------------------------------
        function obj = RunFEAs(obj, config, matProp, loop)
            [~, t2] = size(config.loadingCase);
            for loadcaseIndex = 1:t2
                loadcase = config.loadingCase(loadcaseIndex);
                
                % FE-ANALYSIS
                if (config.w1 ~= 1)
                    u_heat_loadcase  =temperatureFEA_V3(obj, config, matProp,loop,loadcase);
                    obj.U_heatColumn(loadcaseIndex,:)=u_heat_loadcase;
                end
                [UloadCase, obj.maxF, obj.maxU]=FE_elasticV2(obj, config, matProp,loadcase);
                obj.U(loadcaseIndex,:)=UloadCase;
            end
            
            %             obj.c = 0.; % c is the objective. Total strain energy
            %             obj.cCompliance = 0;
            %             obj.cHeat = 0;
        end
        
        % -----------------------------------------
        %
        % Calculate Topology Sensitivity
        %
        % -- FEA data sent to this function.
        %
        % -----------------------------------------
        function obj = CalculateTopologySensitivity(obj, config, matProp, loop)
            elementsInRow = config.nelx+1;
            
            obj.sensitivityHeat(:,:)=0;
            obj.sensitivityElastic(:,:)=0;
            obj.cCompliance=0;
            
            % allow multiple loading cases.
            [~, t2] = size(config.loadingCase);
            
            %             for loadcaseIndex = 1:t2
            % UloadCase= obj.U(loadcaseIndex,:);
            E12_local = 1;
            E33_local = 1;
            
            % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
            count =1;
            for ely = 1:config.nely
                rowMultiplier = ely-1;
                for elx = 1:config.nelx
                    nodes1=[rowMultiplier*elementsInRow+elx;
                        rowMultiplier*elementsInRow+elx+1;
                        (rowMultiplier +1)*elementsInRow+elx+1;
                        (rowMultiplier +1)*elementsInRow+elx];
                    
                    xNodes = nodes1*2-1;
                    yNodes = nodes1*2;
                    NodeNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
                    
                    % if heat objective!!!
                    if (config.w1 ~= 1)
                        U_heat = obj.U_heatColumn(nodes1,:);
                        %                             averageElementTemp = mean2(U_heat); % calculate the average temperature of the 4 nodes
                        KEHeat = matProp.effectiveHeatKEmatrix(  obj.w(ely,elx), config);
                        obj.cHeat =   obj.cHeat           + obj.x(ely,elx)^config.penal*U_heat'*KEHeat*U_heat;
                        
                        % calculate the minim temp sensitivity
                        obj.sensitivityHeat(ely,elx) = -config.penal*obj.x(ely,elx)^(config.penal-1)*U_heat'*KEHeat*U_heat + obj.sensitivityHeat(ely,elx);
                        % obj.g1heat(ely,elx) = obj.x(ely,elx)^(config.penal)*U_heat'*matProp.dKheat*U_heat + obj.g1heat(ely,elx) ;
                    end
                    
                    % Calculate generic KE with 1 as the SIMP density
                     e = elx+(ely-1)*config.nely;
                     if(config.anisotropicMat==1)
                        E12_local=obj.E12(ely,elx);
                        E33_local=obj.E33(ely,elx);
                    end
                    KEgeneric = matProp.getKMatrixTopExxYyyRotVars(config,1,obj.Exx(ely,elx), obj.Eyy(ely,elx),obj.t(ely,elx),obj.w(ely,elx),E12_local, E33_local,e);
                    % Calculate the elastic sensitivity!
                    KEsensitive = KEgeneric*config.penal*obj.x(ely,elx)^(config.penal-1);
                    % Calculate the actual KE
                    KE = KEgeneric*obj.x(ely,elx)^(config.penal);
                    
                    % allow multiple loading cases.
                    Urows = obj.U(:,NodeNumbers);
                    for tt = 1:t2
                        Ue = Urows(tt,:)';
                        obj.sensitivityElastic(ely,elx) =-Ue'*KEsensitive*Ue+ obj.sensitivityElastic(ely,elx);
                        % Calculate the total elastic compliance
                        obj.cCompliance = obj.cCompliance + Ue'*KE*Ue;
                    end
                    
                    count=count+1;
                    
                end % end, loop over x
            end % end, for loop over nely
            %             end % end for loop over load cases.
            
            obj.c=obj.cCompliance*config.w1+obj.cHeat*config.w2;
        end % End Function, CalculateSenstivities
        
        % -----------------------------------------
        %
        % CalculateObjectiveValue
        %
        %  % -- FEA data sent to this function.
        %
        % -----------------------------------------
        function obj = CalculateObjectiveValue(obj, config, matProp, loop,OptimizerLocal)
            elementsInRow = config.nelx+1;
            
            obj.cCompliance=0;
            obj. ExxSysAndSubDiffSummed=0;
            obj. EyySysAndSubDiffSummed=0;
            obj.  ThetaSysAndSubDiffSummed=0;
            
            % allow multiple loading cases.
            [~, t2] = size(config.loadingCase);
            
            %             for loadcaseIndex = 1:t2
            % UloadCase= obj.U(loadcaseIndex,:);
            o = OptimizerLocal;
            if(config.useTargetMesoDensity==1)
                co =   obj. ResponseSurfaceCoefficents;
                sumDensity=0;
            end
            E12_local = 1;
            E33_local = 1;
            
            obj.maxElemStraniEnergy=0;
            % OBJECTIVE FUNCTION
            count =1;
            temp3=obj.Exx*0;
            for ely = 1:config.nely
                rowMultiplier = ely-1;
                for elx = 1:config.nelx
                    nodes1=[rowMultiplier*elementsInRow+elx;
                        rowMultiplier*elementsInRow+elx+1;
                        (rowMultiplier +1)*elementsInRow+elx+1;
                        (rowMultiplier +1)*elementsInRow+elx];
                    
                    xNodes = nodes1*2-1;
                    yNodes = nodes1*2;
                    NodeNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
                    
                    % if heat objective!!!
                    if (config.w1 ~= 1)
                        U_heat = obj.U_heatColumn(nodes1,:);
                        %                             averageElementTemp = mean2(U_heat); % calculate the average temperature of the 4 nodes
                        KEHeat = matProp.effectiveHeatKEmatrix(  obj.w(ely,elx), config);
                        obj.cHeat =   obj.cHeat           + obj.x(ely,elx)^config.penal*U_heat'*KEHeat*U_heat;
                        
                        % calculate the minim temp sensitivity
                        obj.sensitivityHeat(ely,elx) = -config.penal*obj.x(ely,elx)^(config.penal-1)*U_heat'*KEHeat*U_heat + obj.sensitivityHeat(ely,elx);
                        % obj.g1heat(ely,elx) = obj.x(ely,elx)^(config.penal)*U_heat'*matProp.dKheat*U_heat + obj.g1heat(ely,elx) ;
                    end
                    
                    % Calculate generic KE with 1 as the SIMP density
                    e = elx+(ely-1)*config.nely;
                    if(config.anisotropicMat==1)
                        E12_local=obj.E12(ely,elx);
                        E33_local=obj.E33(ely,elx);
                    end
                    KEgeneric = matProp.getKMatrixTopExxYyyRotVars(config,1,obj.Exx(ely,elx), obj.Eyy(ely,elx),obj.t(ely,elx),obj.w(ely,elx), E12_local, E33_local, e);
                    % Calculate the elastic sensitivity!
                    %                     KEsensitive = KEgeneric*config.penal*obj.x(ely,elx)^(config.penal-1);
                    % Calculate the actual KE
                    KE = KEgeneric*obj.x(ely,elx)^(config.penal);
                    
                    % allow multiple loading cases.
                    Urows = obj.U(:,NodeNumbers);
                    elementCompliance=0;
                    for tt = 1:t2
                        Ue = Urows(tt,:)';
                        %                         obj.sensitivityElastic(ely,elx) =-Ue'*KEsensitive*Ue+ obj.sensitivityElastic(ely,elx);
                        % Calculate the total elastic compliance
                        temp =  Ue'*KE*Ue;
                        elementCompliance=elementCompliance+temp;
                        
                    end
                    obj.cCompliance = obj.cCompliance+elementCompliance;
                    
                    % record the largetst element strain energy
                    if( elementCompliance>obj.maxElemStraniEnergy)
                        obj.maxElemStraniEnergy= elementCompliance;
                    end
                    
%                     if(config.useTargetMesoDensity==1)
%                         xxx=obj.Exx(ely,elx)/matProp.E_material1;
%                         yyy=obj.Eyy(ely,elx)/matProp.E_material1;
%                         theta =  obj.t(ely,elx);
%                         
%                         [~, ~,estimateElementDensity] = o.CalculateDensitySensitivityandRho(xxx,yyy,theta,obj.ResponseSurfaceCoefficents,config,matProp);
%                         
%                         estimateElementDensity= min(max(estimateElementDensity,0.05),1);%1 is max, 0.05 is min
%                         eleDensity = obj.x(ely,elx)*estimateElementDensity;
%                         sumDensity =sumDensity+eleDensity;
%                         temp3(ely,elx) = eleDensity;
%                     end
                    
                    
                    count=count+1;
                    
                end % end, loop over x
            end % end, for loop over nely
            %             end % end for loop over load cases.
            
            % Calculate how well the consistency constraints are working.
            if(config.macro_meso_iteration>1)
                obj.ExxSysAndSubDiffSummed=sum(sum(obj.x.*((obj.Exx-obj.ExxSub).^2).^(1/2))); % make sure it is absolute value, by square, then sqrt
                obj.EyySysAndSubDiffSummed=sum(sum(obj.x.*((obj.Eyy-obj.EyySub).^2).^(1/2)));% sum(sum(sqrt(obj.x.*(obj.Eyy-obj.EyySub).^2)));
                
                % if Exx = Eyy, then don't count in the  theta calcualtions
                logicArray = abs(obj.Eyy-obj.Exx)./obj.Exx<=0.02;
                logicArray=1-logicArray;
                obj.ThetaSysAndSubDiffSummed=sum(sum(logicArray.*obj.x.*((obj.t-obj.thetaSub).^2).^(1/2)));%sum(sum(sqrt(obj.x.*(obj.t-obj.thetaSub).^2)));
            else
                obj.ExxSysAndSubDiffSummed=0;
                obj.EyySysAndSubDiffSummed  =0;
                obj.ThetaSysAndSubDiffSummed=0;
            end
          
            
            if(config.useTargetMesoDensity==1)
                 theta = obj.t;
%                logic2 = theta>pi/4;
%                logic3 = 0<theta<pi/4;
%                  logic4 = theta<0;
%                theta(logic2)=theta(logic2)-pi/4;
%                theta(logic3)=pi/4-theta(logic3);
%                 theta(logic4)=-pi/4-theta(logic4);
                  
                [~, ~,rhoValue] = o.CalculateDensitySensitivityandRho(obj.Exx/matProp.E_material1,obj.Eyy/matProp.E_material1,theta,obj.x ,obj.ResponseSurfaceCoefficents,config,matProp,obj.densityOffsetArray);
              
                temp2 = sum(sum(rhoValue));
                sumDensity=temp2/(config.nelx*config.nely*config.totalVolume);
%                 sumDensity = sumDensity/(config.nelx*config.nely*config.totalVolume);
            else
                sumDensity=0;
            end
            
            obj.averageMesoDensity=sumDensity;
            %                obj.averageMesoDensity=mesoDensity;
            obj.c=obj.cCompliance*config.w1+obj.cHeat*config.w2;
            if(obj.c<0)
                t = 'something really wrong'
            end
        end
        
        function obj = AddDataToStoreOptimizationVarArray(obj,config)
            obj.storeOptimizationVar = [obj.storeOptimizationVar;...
                obj.c, ...
                obj.cCompliance,...
                obj.cHeat, ...
                obj.currentVol1Fraction, ...
                obj.currentVol2Fraction, ...
                sum(sum(obj.x)),...
                config.targetAvgExxEyy, ...
                obj.actualAverageE, ...
                obj.averageMesoDensity...
                obj.ExxSysAndSubDiffSummed...
                obj.EyySysAndSubDiffSummed...
                obj.ThetaSysAndSubDiffSummed...
                ];
        end
        
        
        % -----------------------------------------
        %
        % CalculateTMaterialGradientSensitivity
        %
        %  % -- FEA data sent to this function.
        %
        % -----------------------------------------
        function obj = CalculateMaterialGradientSensitivity(obj, config, matProp, loop)
            elementsInRow = config.nelx+1;
            
            % allow multiple loading cases.
            [~, t2] = size(config.loadingCase);
            
            obj.sensitivityHeat(:,:)=0;
            obj.sensitivityElastic(:,:)=0;
            
            
            %             for loadcaseIndex = 1:t2
            %                 UloadCase=obj.U(loadcaseIndex,:);
            %                 count =1;
            
            for ely = 1:config.nely
                rowMultiplier = ely-1;
                for elx = 1:config.nelx
                    nodes1=[rowMultiplier*elementsInRow+elx;
                        rowMultiplier*elementsInRow+elx+1;
                        (rowMultiplier +1)*elementsInRow+elx+1;
                        (rowMultiplier +1)*elementsInRow+elx];
                    
                    xNodes = nodes1*2-1;
                    yNodes = nodes1*2;
                    NodeNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
                    
                    %                         Ue = UloadCase(NodeNumbers)';
                    
                    % Calculate the elastic sensitivity!
                    KEsensitive = matProp.getKMatrixGradientMaterialSensitivity(config,obj.x(ely,elx),obj.Exx(ely,elx), obj.Eyy(ely,elx),obj.t(ely,elx));
                    
                    % allow multiple loading cases.
                    Urows = obj.U(:,NodeNumbers);
                    for tt = 1:t2
                        Ue = Urows(tt,:)';
                        obj.sensitivityElastic(ely,elx) =Ue'*KEsensitive*Ue+ obj.sensitivityElastic(ely,elx);
                        % Calculate the total elastic compliance
                        %                         obj.cCompliance = obj.cCompliance + Ue'*KE*Ue;
                    end
                    
                    
                    %                         obj.sensitivityElastic(ely,elx) =Ue'*KEsensitive*Ue+ obj.sensitivityElastic(ely,elx);
                    
                    
                    % if heat objective!!!
                    if (config.w1 ~= 1)
                        U_heat = obj.U_heatColumn(nodes1,:);
                        %                             averageElementTemp = mean2(U_heat); % calculate the average temperature of the 4 nodes
                        KEHeat = matProp.effectiveHeatKEmatrix(  obj.w(ely,elx), config);
                        obj.cHeat =   obj.cHeat  + obj.x(ely,elx)^config.penal*U_heat'*KEHeat*U_heat;
                        
                        % calculate the minim temp sensitivity
                        %                             obj.temp2(ely,elx) = -config.penal*obj.x(ely,elx)^(config.penal-1)*U_heat'*KEHeat*U_heat + obj.temp2(ely,elx);
                        obj.sensitivityHeat(ely,elx) = obj.x(ely,elx)^(config.penal)*U_heat'*matProp.dKheat*U_heat + obj.sensitivityHeat(ely,elx) ;
                    end
                    
                    %                         count=count+1;
                end % end, loop over x
            end % end, for loop over y
            
            %             end % end for loop over load cases.
        end % End Function, CalculateMaterialGradientSensitivity
        
        
        
        
        % -----------------------------------------
        %
        % CalculateExxEyySensitivity
        %
        %  % -- FEA data sent to this function.
        %
        % -----------------------------------------
        function obj = CalculateExxEyySensitivity(obj, config, matProp, loop)
            elementsInRow = config.nelx+1;
            obj.sensitivityHeat(:,:)=0;
            obj.sensitivityElastic(:,:)=0;
            obj.sensitivityElasticPart2(:,:) =0;
            
            % allow multiple loading cases.
            [~, t2] = size(config.loadingCase);
            
            %             for loadcaseIndex = 1:t2
            %                 UloadCase= obj.U(loadcaseIndex,:);
            E12_local = 1;
            E33_local = 1;
            %  SENSITIVITY ANALYSIS
            for y = 1:config.nely
                rowMultiplier = y-1;
                for xx = 1:config.nelx
                    nodes1=[rowMultiplier*elementsInRow+xx;
                        rowMultiplier*elementsInRow+xx+1;
                        (rowMultiplier +1)*elementsInRow+xx+1;
                        (rowMultiplier +1)*elementsInRow+xx];
                    
                    xNodes = nodes1*2-1;
                    yNodes = nodes1*2;
                    NodeNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
                    
                    URows = obj.U(:,NodeNumbers);
                    
                    % Get the sensitivity K matrix
                    % Set Exx = 1, Eyy = 0 to get sensitivity
                    if(config.useRinOrthMaterialModel==0)                  
                        KExx = matProp.getKMatrixTopExxYyyRotVars(config,obj.x(y,xx),1, 0,obj.t(y,xx),E12_local, E33_local,[]);                    
                        % Set Eyy= 1,  Eyy = 0,to get sensitivity
                        KEyy = matProp.getKMatrixTopExxYyyRotVars(config,obj.x(y,xx),0, 1,obj.t(y,xx),E12_local, E33_local,[]);
                    else
                        topDensity=obj.x(y,xx);
                        Exx_local=obj.Exx(y,xx);
                        Eyy_local=obj.Eyy(y,xx);
                        rotation=obj.t(y,xx);
                        
                        mode =1;
                        KExx = matProp.getKMatrixSensitivityTopExxYyyRotVarsWithRshear(config,topDensity,Exx_local, Eyy_local,rotation,mode)  ;
                        mode =2;
                        KEyy = matProp.getKMatrixSensitivityTopExxYyyRotVarsWithRshear(config,topDensity,Exx_local, Eyy_local,rotation,mode);
                    end
                    
                    
                    % allow multiple loading cases.
                    tempSensi = 0;
                    tempSensiPart2=0;
                    for i = 1:t2
                        Ucase = URows(i,:)';
                        tempSensi= tempSensi+Ucase'*KExx*Ucase;
                        tempSensiPart2 = tempSensiPart2++Ucase'*KEyy*Ucase;
                    end
                    
                    obj.sensitivityElastic(y,xx) = tempSensi + obj.sensitivityElastic(y,xx);
                    obj.sensitivityElasticPart2(y,xx) = tempSensiPart2 + obj.sensitivityElasticPart2(y,xx);
                    
                end % end, loop over x
            end % end, for loop over nely
            %             end % end for loop over load cases.
            
        end % End Function, CalculateOthogonalDistributionSensitivity
        
        % -----------------------------------------
        %
        % CalculateANISOTROPICSensitivity
        %
        %  % -- FEA data sent to this function.
        %
        % -----------------------------------------
        function obj = CalculateANISOTROPICSensitivity(obj, config, matProp, loop)
            elementsInRow = config.nelx+1;
            obj.sensitivityHeat(:,:)=0;
            obj.sensitivityElastic(:,:)=0;
            obj.sensitivityElasticPart2(:,:) =0;
            obj. sensitivityElasticE12(:,:)=0;
            obj. sensitivityElasticE33(:,:)=0;
            
            % allow multiple loading cases.
            [~, t2] = size(config.loadingCase);
            
            %             for loadcaseIndex = 1:t2
            %                 UloadCase= obj.U(loadcaseIndex,:);
%             config,topDensity,Exx, Eyy,rotation,material1Fraction,E12, E33,e)
            thetaTemp = 1;
            wTemp = 1;
            simpTemp =1;
            eTemp=1;
             KE_xx = matProp.getKMatrixTopExxYyyRotVars(config,simpTemp,1, 0,thetaTemp, wTemp,0,0,eTemp);
             KE_yy = matProp.getKMatrixTopExxYyyRotVars(config,simpTemp,0, 1,thetaTemp, wTemp,0,0,eTemp);
             KE_e12 = matProp.getKMatrixTopExxYyyRotVars(config,simpTemp,0, 0,thetaTemp, wTemp,1,0,eTemp);
             KE_e33 = matProp.getKMatrixTopExxYyyRotVars(config,simpTemp,0, 0,thetaTemp, wTemp,0,1,eTemp);
              
              
            %  SENSITIVITY ANALYSIS
            for y = 1:config.nely
                rowMultiplier = y-1;
                for xx = 1:config.nelx
                    nodes1=[rowMultiplier*elementsInRow+xx;
                        rowMultiplier*elementsInRow+xx+1;
                        (rowMultiplier +1)*elementsInRow+xx+1;
                        (rowMultiplier +1)*elementsInRow+xx];
                    
                    xNodes = nodes1*2-1;
                    yNodes = nodes1*2;
                    NodeNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
                    
                    URows = obj.U(:,NodeNumbers);
                    
                    % Get the sensitivity K matrix
                    % Set Exx = 1, Eyy = 0 to get sensitivity
                    p=obj.x(y,xx)^config.penal;
                    
                    KE_xx_temp = KE_xx*p;
                    KE_yy_temp = KE_yy*p;
                    KE_e12_temp = KE_e12*p;
                    KE_e33_temp = KE_e33*p;
                    
                    
                    
                    % allow multiple loading cases.
                    tempSensi = 0;
                    tempSensiPart2=0;
                    tempSensi_E12=0;
                    tempSensi_E33=0;
                    for i = 1:t2
                        Ucase = URows(i,:)';
                        tempSensi= tempSensi+Ucase'*KE_xx_temp*Ucase;
                        tempSensiPart2 = tempSensiPart2+Ucase'*KE_yy_temp*Ucase;
                        tempSensi_E12 = tempSensi_E12+Ucase'*KE_e12_temp*Ucase;
                        tempSensi_E33=tempSensi_E33+Ucase'*KE_e33_temp*Ucase;
                    end
                    
                    obj.sensitivityElastic(y,xx) = tempSensi;% + obj.sensitivityElastic(y,xx);
                    obj.sensitivityElasticPart2(y,xx) = tempSensiPart2;% + obj.sensitivityElasticPart2(y,xx);
                    obj. sensitivityElasticE12(y,xx)=tempSensi_E12;
                    obj. sensitivityElasticE33(y,xx)=tempSensi_E33;
                    
                end % end, loop over x
            end % end, for loop over nely
            %             end % end for loop over load cases.
            
        end % End Function, CalculateOthogonalDistributionSensitivity
        
        % ------------------------------------------------------------
        %
        % CalculateSensitiviesMesoStructure no periodic
        %
        % ------------------------------------------------------------
%         function obj = CalculateSensitiviesMesoStructureNoPeriodic(obj, config, matProp, loop,macroElemProps, U)
%             
%             doplot = 0;
%             doplotfinal = 0;
%             if(doplot ==1 || doplotfinal==1)
%                 
%                 p = plotResults;
%                 figure(1);
%             end
%             
%             % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
%             obj.c = 0.; % c is the objective. Total strain energy
%             obj.cCompliance = 0;
%             obj.cHeat = 0;
%             
%             [~, t2] = size(config.loadingCase);
%             % allow multiple loading cases.
%             
%             
%             
%             for loadcaseIndex = 1:t2
%                 Ucase = U(:,loadcaseIndex);
%                 loadcase = config.loadingCase(loadcaseIndex);
%                 ne = config.nelx*config.nely; % number of elements
%                 for e = 1:ne
%                     
%                     % loop over local node numbers to get their node global node numbers
%                     nodes1 = obj.IEN(e,:);
%                     [elx,ely]= obj.GivenNodeNumberGetXY(e);
%                     
%                     xNodes = nodes1*2-1;
%                     yNodes = nodes1*2;
%                     dofNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
%                     
%                     Ue = Ucase(dofNumbers);
%                     
%                     % U_heat = obj.U_heatColumn(nodes1,:);
%                     %averageElementTemp = mean2(U_heat); % calculate the average temperature of the 4 nodes
%                     
%                     % Get the element K matrix for this partiular element
%                     KE = matProp.effectiveElasticKEmatrix(  obj.w(ely,elx),config,[]);
%                     
%                     % KEHeat = matProp.effectiveHeatKEmatrix(  obj.w(ely,elx), config);
%                     % Dmaterial = matProp.calculateEffectiveConstitutiveEquation( obj.w(ely,elx), config);
%                     %                 config.nelx
%                     % Find the elastic strain
%                     elasticStrain = obj.B*Ue;
%                     
%                     % term1 = transpose(Ue)*KE*Ue*obj.x(ely,elx)^(config.penal-1)*config.penal;
%                     
%                     term1 = transpose(Ue)*KE*Ue*obj.x(ely,elx)^(config.penal-1)*config.penal;
%                     
%                     %                      term1_method2 = (eye(3)-elasticStrain);
%                     %  term2 = 0;
%                     % term3= 0;
%                     
%                     % Sum the elastic compliance terms.
%                     % total = (term1 + term2 + term3);
%                     obj.temp1(ely,elx) = term1+obj.temp1(ely,elx);
%                     
%                     if(doplot ==1)
%                         %                     if(mod(e,10) ==0)
%                         
%                         p.PlotArrayGeneric(obj.temp1, 'plotting sensitivities while running. ')
%                         drawnow
%                         %                     end
%                     end
%                     % calculate the minim temp sensitivity
%                     % obj.temp2(ely,elx) = -config.penal*obj.x(ely,elx)^(config.penal-1)*U_heat'*KEHeat*U_heat;
%                 end
%                 
%                 % Do final plot
%                 if(doplotfinal ==1)
%                     subplot(2,2,3);
%                     p.PlotArrayGeneric(obj.temp1, 'final plotting sensitivities after running. ')
%                     drawnow
%                 end
%                 %             end
%             end % end loading cases
%             
%             
%             %                obj.temp1(ely,elx)  =    obj.temp1(ely,elx) /t2; % average the cases
%             obj.temp1  =    obj.temp1 /t2; % average the cases
%         end % end CalculateSensitiviesMesoStructureNoPeriodic
        
        
        
        % ------------------------------------------------------------
        % GetHomogenizedProperties
        % ------------------------------------------------------------
        %
        % For meso strudcture design HOmogenize the matrix and get K_h
        % Return the K_matrix that discribes this meso structure's
        % properties
        % ------------------------------------------------------------
        %         function macroElemProps = GetHomogenizedProperties(obj, config,homgSettings, matProp, loop,macroElemProps)
        %
        %             [macroElemProps.D_homog]=FE_elasticV2_homgonization(obj, config, matProp);
        %
        %         end % end CalculateSensitiviesMesoStructure
        
        
        
        
        %%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [dcn]=check(obj, nelx,nely,rmin,x,dc)
%             dcn=zeros(nely,nelx);
%             for i = 1:nelx
%                 for j = 1:nely
%                     sum=0.0;
%                     for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
%                         for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
%                             fac = rmin-sqrt((i-k)^2+(j-l)^2);
%                             sum = sum+max(0,fac);
%                             dcn(j,i) = dcn(j,i) + max(0,fac)*x(l,k)*dc(l,k);
%                         end
%                     end
%                     dcn(j,i) = dcn(j,i)/(x(j,i)*sum);
%                 end
%             end

dcn=imgaussfilt(dc,rmin );
        end
        
        %  -----------------------------------
        % Flip Exx and Eyy so that theta is always positive
        % %  -----------------------------------
        function [obj]=FlipOrientation(obj,config)
            
            for i = 1:config.nelx
                for j = 1:config.nely
                       localTheta=obj.t(j,i);
                    if(localTheta<0)
                     
                        Exxtemp1=obj.Exx(j,i);
                        EyyTempSee= obj.Eyy(j,i);
                        ThetaTempSee=obj.t(j,i);
                        obj.Exx(j,i)=obj.Eyy(j,i);
                        obj.Eyy(j,i)=Exxtemp1;
                        obj.t(j,i)=obj.t(j,i)+pi/2;
                        
                        % Also flip the meso level values or else the
                        % optimization will be crazy!
                        ExxSubtemp2 = obj.ExxSub(j,i);
                         EyysubTempSee= obj.EyySub(j,i);
                        ThetasubTempSee=obj.thetaSub(j,i);
                        obj.ExxSub(j,i)=obj.EyySub(j,i);   % Sub system copies of design var
                        obj.EyySub(j,i)=ExxSubtemp2;
%                         obj.thetaSub(j,i)= min(obj.thetaSub(j,i)+pi/2,pi/2);
                           obj.thetaSub(j,i)= obj.thetaSub(j,i)+pi/2;
                        
                        
                    end
                    
                    if(localTheta>pi/2)
                     
                        Exxtemp1=obj.Exx(j,i);
                        EyyTempSee= obj.Eyy(j,i);
                        ThetaTempSee=obj.t(j,i);
                        obj.Exx(j,i)=obj.Eyy(j,i);
                        obj.Eyy(j,i)=Exxtemp1;
                        obj.t(j,i)=obj.t(j,i)-pi/2;
                        
                        % Also flip the meso level values or else the
                        % optimization will be crazy!
                        ExxSubtemp2 = obj.ExxSub(j,i);
                         EyysubTempSee= obj.EyySub(j,i);
                        ThetasubTempSee=obj.thetaSub(j,i);
                        obj.ExxSub(j,i)=obj.EyySub(j,i);   % Sub system copies of design var
                        obj.EyySub(j,i)=ExxSubtemp2;
%                         obj.thetaSub(j,i)= min(obj.thetaSub(j,i)-pi/2,0);
                          obj.thetaSub(j,i)= obj.thetaSub(j,i)-pi/2;
                        
                        
                    end
                end
            end
        end
        
        
        
        % -----------------------------
        % GenerateStartingMesoDesign
        %
        % Geneate the x (density) values for the Meso design problem.
        % Several methods exist.
        % -----------------------------
        function [obj]= GenerateStartingMesoDesign(obj,mesoConfig,macroElementProperties)
            method =mesoConfig.mesoDesignInitalConditions;
            
            if(method ==1)
                obj.x(1:mesoConfig.nely,1:mesoConfig.nelx) = mesoConfig.totalVolume; % artificial density of the elements
                % method 1, randome values. Does not seem to be working well.
                obj.x(1:mesoConfig.nely,1:mesoConfig.nelx) = randi([0, round(mesoConfig.totalVolume*100)],mesoConfig.nely,mesoConfig.nelx)/100; % artificial density of the elements, can not be unifrom or else sensitivity will be 0 everywhere.
            elseif(method ==2)
                % method 2, box of empty in the middle.
                obj.x(1:mesoConfig.nely,1:mesoConfig.nelx) = ones(mesoConfig.nely,mesoConfig.nelx);
                midY = round(mesoConfig.nely/2);
                midX = round(mesoConfig.nelx/2);
                ratio = mesoConfig.nelx/mesoConfig.nely;
                vEmpty = mesoConfig.nelx*mesoConfig.nely-mesoConfig.totalVolume*mesoConfig.nelx*mesoConfig.nely;
                dimY = floor(sqrt(vEmpty/ratio));
                yStart = midY-floor(dimY/2);
                dimX =  floor(ratio*dimY);
                xStart = midX-floor(dimX/2);
                obj.x(yStart:yStart+dimY-1,xStart:xStart+dimX-1)= ones(dimY,dimX)*0.01;
                
                
            elseif(method ==3)
                % method 3, circle in the moddle
                obj.x(1:mesoConfig.nely,1:mesoConfig.nelx) = ones(mesoConfig.nely,mesoConfig.nelx);
                %            obj.x(1:mesoConfig.nely,1:mesoConfig.nelx) = randi([1, round(mesoConfig.totalVolume*100)],mesoConfig.nely,mesoConfig.nelx)/100; % artificial density of the elements, can not be unifrom or else sensitivity will be 0 everywhere.
                
                midY = round(mesoConfig.nely/2);
                midX = round(mesoConfig.nelx/2);
%                 mesoConfig.totalVolume=0.5;%*mesoConfig.nely*mesoConfig.nelx;
                radius = sqrt((mesoConfig.totalVolume*mesoConfig.nelx*mesoConfig.nely-mesoConfig.nelx*mesoConfig.nely)/(-pi));
                for i = 1:mesoConfig.nelx
                    for j = 1:mesoConfig.nely
                        %                 if sqrt((i-mesoConfig.nelx/2-0.5)^2+(j-mesoConfig.nely/2-0.5)*2) < min(mesoConfig.nelx,mesoConfig.nely)/3
                        %                     obj.x(j,i) = mesoConfig.totalVolume/2;
                        %                 end
                        d = sqrt((i-midX)^2+(j-midY)^2);
                        if(d<radius)
                            obj.x(j,i)= 0.01;
                        end
                    end
                end
                
            elseif(method ==4)
                % method 4, many circle holes
                obj.x(1:mesoConfig.nely,1:mesoConfig.nelx) = ones(mesoConfig.nely,mesoConfig.nelx);
                numHolesX = 7;
                numHolesY =7;
                totalHoles =numHolesX*numHolesY;
                XholeCenters = 1:mesoConfig.nelx/(numHolesX+1):mesoConfig.nelx;
                YholeCenters = 1:mesoConfig.nely/(numHolesY+1):mesoConfig.nely;
                
                radius = sqrt((mesoConfig.totalVolume*mesoConfig.nelx*mesoConfig.nely-mesoConfig.nelx*mesoConfig.nely)/(-pi*totalHoles));
                for i = 1:mesoConfig.nelx
                    for j = 1:mesoConfig.nely
                        %                 if sqrt((i-mesoConfig.nelx/2-0.5)^2+(j-mesoConfig.nely/2-0.5)*2) < min(mesoConfig.nelx,mesoConfig.nely)/3
                        %                     obj.x(j,i) = mesoConfig.totalVolume/2;
                        %                 end
                        for mm = 1:(1+numHolesX)
                            for nn = 1:(1+numHolesY)
                                midX= XholeCenters(mm);
                                midY= YholeCenters(nn);
                                
                                d = sqrt((i-midX)^2+(j-midY)^2);
                                if(d<radius)
                                    obj.x(j,i)=  0.01;
                                end
                            end
                        end
                    end
                end
            elseif(method ==5)
                % method 5, small hole in middle
                obj.x(1:mesoConfig.nely,1:mesoConfig.nelx) = ones(mesoConfig.nely,mesoConfig.nelx);
                midY = round(mesoConfig.nely/2);
                midX = round(mesoConfig.nelx/2);
                radius = 3;
                for i = 1:mesoConfig.nelx
                    for j = 1:mesoConfig.nely
                        %                 if sqrt((i-mesoConfig.nelx/2-0.5)^2+(j-mesoConfig.nely/2-0.5)*2) < min(mesoConfig.nelx,mesoConfig.nely)/3
                        %                     obj.x(j,i) = mesoConfig.totalVolume/2;
                        %                 end
                        d = sqrt((i-midX)^2+(j-midY)^2);
                        if(d<radius)
                            obj.x(j,i)=  0.01;
                        end
                        
                    end
                end
                
                % Randome circles
            elseif(method ==6)
                % method 4, many randome circle holes
                numHolesX =6;
                numHolesY =6;
                totalHoles =numHolesX*numHolesY;
                XholeCenters = randi([0, mesoConfig.nelx],1,numHolesX+1);
                YholeCenters =  randi([0, mesoConfig.nely],1,numHolesY+1);
                
                obj.x(1:mesoConfig.nely,1:mesoConfig.nelx) = ones(mesoConfig.nely,mesoConfig.nelx);
                
                radius = sqrt((mesoConfig.totalVolume*mesoConfig.nelx*mesoConfig.nely-mesoConfig.nelx*mesoConfig.nely)/(-pi*totalHoles));
                for i = 1:mesoConfig.nelx
                    for j = 1:mesoConfig.nely
                        %                 if sqrt((i-mesoConfig.nelx/2-0.5)^2+(j-mesoConfig.nely/2-0.5)*2) < min(mesoConfig.nelx,mesoConfig.nely)/3
                        %                     obj.x(j,i) = mesoConfig.totalVolume/2;
                        %                 end
                        for mm = 1:(1+numHolesX)
                            for nn = 1:(1+numHolesY)
                                midX= XholeCenters(mm);
                                midY= YholeCenters(nn);
                                
                                d = sqrt((i-midX)^2+(j-midY)^2);
                                if(d<radius)
                                    obj.x(j,i)=  0.01;
                                end
                            end
                            
                        end
                    end
                end
            end
        end
        
        % -----------------------------
        % GenerateMesoStructureCoordinationMask
        %
        % Geneate the 1 or 0 mask in order to coordinate the meso structure
        % designs so that they actually touch.
        % -----------------------------
        function [obj]= GenerateMesoStructureCoordinationMask(obj,mesoConfig,macroElementProps)
            config = mesoConfig;
              obj.mesoStructNTCmask = zeros(config.nely,config.nelx);
            if (config.coordinateMesoBoundaries==1 && config.macro_meso_iteration>1)
                   mm_iteration = config.macro_meso_iteration;
                   
                   % Get the density field
                   outname = sprintf('./out%i/SIMPdensityfield%i.csv',config.iterationNum,mm_iteration);
                   xx = csvread(outname);              
                % preform the logic tests  to see if we need to add material     
                add_top=0;
                add_bottom=0;
                add_left=0;
                add_right=0;
                
                xCurrent=macroElementProps.xPos;
                xRight=xCurrent+1;
                xLeft=xCurrent-1;
                
                yCurrent=macroElementProps.yPos;
                yUp=yCurrent+1;
                yDown=yCurrent-1;    
                
                % -----------------------------
                % Check Right
                % -------------------
                if(xRight<=config.nelxMacro)
                    
                    density = xx(yCurrent, xRight);
                    if(density>config.voidMaterialDensityCutOff)
                        add_right=1;
                    end
                end
                
                
                % -----------------------------
                % Check Left
                % -------------------
                if(xLeft>0)
                    
                    density = xx(yCurrent, xLeft);
                    if(density>config.voidMaterialDensityCutOff)
                        add_left=1;
                    end
                end
                
                
                % -----------------------------
                % Check down
                % -------------------
                if(yDown>0)
                    
                    density = xx(yDown, xCurrent);
                    if(density>config.voidMaterialDensityCutOff)
                        add_bottom=1;
                    else
%                         % Check Right, then Check left; checking for empty
%                         if(xRight<config.nelxMacro)
%                             densityxRight= xx(yCurrent, xRight);
%                             if(densityxRight<config.voidMaterialDensityCutOff)
%                                 cut_bottomRight=1;
%                             end
%                         end
%                         
%                         if(xLeft>0)
%                             densityxLeft= xx(yCurrent, xLeft);
%                             if(densityxLeft<config.voidMaterialDensityCutOff)
%                                 cut_bottomLeft=1;
%                             end
%                         end
                    end
                end
                
                
                % -----------------------------
                % Check up
                % -------------------
                if(yUp<=config.nelyMacro)
                    
                    densityUp = xx(yUp, xCurrent);
                    if(densityUp>config.voidMaterialDensityCutOff)
                        add_top=1;
                    else
%                         % Check Right, then Check left; checking for empty
%                         if(xRight<config.nelxMacro)
%                             densityxRight= xx(yCurrent, xRight);
%                             if(densityxRight<config.voidMaterialDensityCutOff)
%                                 cut_topRight=1;
%                             end
%                         end
%                         
%                         if(xLeft>0)
%                             densityxLeft= xx(yCurrent, xLeft);
%                             if(densityxLeft<config.voidMaterialDensityCutOff)
%                                 cut_topLeft=1;
%                             end
%                         end
                    end
                end
                
                % --------------------------
                % Add the boundaries to the mask
                % --------------------------
                
                
                
                % GET the saved element to XY position map (needed for x and w vars retrival)
%                 outname = sprintf('./out%i/elementXYposition%i.csv',config.iterationNum,config.macro_meso_iteration);
%                 elementXYpositionMacro=csvread(outname);
%                 results = elementXYpositionMacro(macroEleProps.elementNumber,:);
%                 macroEleProps.yPos = results(1);
%                 macroEleProps.xPos = results(2);
                settings=config;
                settings.macro_meso_iteration =config.macro_meso_iteration-1;
                settings.nelx = settings.nelxMacro;
                settings.nely = settings.nelyMacro;

                numrows=1;
                sumIndex = xCurrent+yCurrent;

                doPlot =1;
                if(doPlot ==1)
                    p = plotResults;
                end

                % switch between even and odd ever macro, meso iteration
                evenOddTarget = mod(config.macro_meso_iteration,2);
                % pick ever other one.
                if(mod(sumIndex,2)==evenOddTarget)


                    % Get the design above, and get the previous
                    % iteration's design


                    if (add_top==1)
                        elementNumber =obj. GetNodeNumberGivenXY(settings, xCurrent,yUp);
                        [xAdjacent]= GetMesoUnitCellDesignFromCSV(settings,elementNumber);
                        obj.mesoStructNTCmask(end-numrows+1:end,:)=xAdjacent(end-numrows+1:end,:);

                        if(doPlot ==1)
                            subplot(3,3,2);
                            p.PlotArrayGeneric(xAdjacent,'up');
                        end
                    end


                    if (add_bottom==1)
                        elementNumber = obj.GetNodeNumberGivenXY(settings, xCurrent,yDown);
                        [xAdjacent]= GetMesoUnitCellDesignFromCSV(settings,elementNumber);
                        obj.mesoStructNTCmask(1:numrows,:)=xAdjacent(1:numrows,:);

                        if(doPlot ==1)
                            subplot(3,3,8);
                            p.PlotArrayGeneric(xAdjacent,'down');
                        end
                    end

                    if (add_right==1)
                        elementNumber =obj. GetNodeNumberGivenXY(settings, xRight,yCurrent);
                        [xAdjacent]= GetMesoUnitCellDesignFromCSV(settings,elementNumber);
                        obj.mesoStructNTCmask(:,end-numrows+1:end)=xAdjacent(:,end-numrows+1:end);

                        if(doPlot ==1)
                            subplot(3,3,6);
                            p.PlotArrayGeneric(xAdjacent,'right');
                        end
                    end

                    if(add_left==1)
                        elementNumber =obj. GetNodeNumberGivenXY(settings, xLeft,yCurrent);
                        [xAdjacent]= GetMesoUnitCellDesignFromCSV(settings,elementNumber);
                        obj.mesoStructNTCmask(:,1:numrows)=xAdjacent(:,1:numrows);

                        if(doPlot ==1)
                            subplot(3,3,4);
                            p.PlotArrayGeneric(xAdjacent,'left');
                        end
                    end
                    
                    %obj.mesoStructNTCmask( obj.mesoStructNTCmask>config.voidMaterialDensityCutOff)=1;
                     %obj.mesoStructNTCmask( obj.mesoStructNTCmask<config.voidMaterialDensityCutOff)=0;
                    if(doPlot ==1)
                        subplot(3,3,5);
                        p.PlotArrayGeneric(obj.mesoStructNTCmask,'mask');
                    end
                else
                    %------------------------------
                    % Add boundaries from previous self design
                    %------------------------------
                    elementNumber =obj. GetNodeNumberGivenXY(settings, xCurrent,yCurrent);
                    [xAdjacent]= GetMesoUnitCellDesignFromCSV(settings,elementNumber);

                    % All the 4 sides combined
                    obj.mesoStructNTCmask(end-numrows+1:end,:)=xAdjacent(end-numrows+1:end,:);
                    obj.mesoStructNTCmask(1:numrows,:)=xAdjacent(1:numrows,:);
                    obj.mesoStructNTCmask(:,end-numrows+1:end)=xAdjacent(:,end-numrows+1:end);
                    obj.mesoStructNTCmask(:,1:numrows)=xAdjacent(:,1:numrows);
                    
                      %obj.mesoStructNTCmask( obj.mesoStructNTCmask>config.voidMaterialDensityCutOff)=1;
                %obj.mesoStructNTCmask( obj.mesoStructNTCmask<config.voidMaterialDensityCutOff)=0;

                    if(doPlot ==1)
                        subplot(2,1,1);
                        p.PlotArrayGeneric(obj.mesoStructNTCmask,'mask');

                        subplot(2,1,2);
                        p.PlotArrayGeneric(xAdjacent,'old');
                    end

                end
              
                
            end
        end
        
        % -----------------------------
        % AddCoordinationMaskToSensitivies
        %
        % -----------------------------
        function [newSensitivities]= AddCoordinationMaskToSensitivies(obj,mesoConfig,macroElementProperties)
            config = mesoConfig;
            if (config.coordinateMesoBoundaries==1 && config.macro_meso_iteration>1)
%                  offset =median(median(obj.dc));
                   offset =min(min(obj.dc));
                   offset=offset*0.5
                 offset=offset*(config.macro_meso_iteration-1); % INcrease the strength each time
%                 offset =min(min(obj.dc));
               newSensitivities=obj.dc+ obj.mesoStructNTCmask*offset;
            else
                newSensitivities=obj.dc;
            end
        end
        
    end % End Methods
end % End Class DesignVars