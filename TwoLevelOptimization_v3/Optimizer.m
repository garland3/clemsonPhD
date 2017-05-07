classdef Optimizer
    % Optimizer - class contains only the code that actually optimizes.
    %
    
    properties
    end
    
    methods
        % -------------------------------
        % TOPOOLOGY, SIMP METHOD
        % -------------------------------
        function DV = OptimizeTopology(obj,DV, config, matProp,masterloop)
            DV = DV.CalculateTopologySensitivity(config, matProp, masterloop);
            % normalize the sensitivies  by dividing by their max values.
            if (config.w1 ~= 1) % if we are using the heat objective
                temp1Max =-1* min(min(DV.sensitivityElastic));
                DV.sensitivityElastic = DV.sensitivityElastic/temp1Max;
                temp2Max = -1* min(min(DV.sensitivityHeat));
                DV.sensitivityHeat = DV.sensitivityHeat/temp2Max;
                DV.dc = config.w1*DV.sensitivityElastic+config.w2*DV.sensitivityHeat; % add the two sensitivies together using their weights
            else
                DV.dc = config.w1*DV.sensitivityElastic;
            end
            % FILTERING OF SENSITIVITIES
            [DV.dc]   = DV.check( config.nelx, config.nely,config.rmin,DV.x,DV.dc);
            % DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
            moveLimit=0.1;
            [DV.x]    = OC( config.nelx, config.nely,DV.x,config.totalVolume,DV.dc, DV, config,moveLimit);
        end
        
        % ----------------------------------
        % VOLUME FRACTION OPTIMIZATION
        % ----------------------------------
        function DV = OptimizeVolumeFraction(obj,DV,config, matProp, masterloop)
            DV = DV.CalculateMaterialGradientSensitivity(config, matProp, masterloop);
            
            DV =  DV.CalculateVolumeFractions(config,matProp);
            
            totalVolLocal = DV.currentVol1Fraction+ DV.currentVol2Fraction;
            fractionCurrent_V1Local = DV.currentVol1Fraction/totalVolLocal;
            targetFraction_v1 = config.v1/(config.v1+config.v2);
            
            % Normalize the sensitives.
            if (config.w1 ~= 1) % if we are using the heat objective
                temp1Max = max(max(abs(DV.sensitivityElastic)));
                DV.sensitivityElastic = DV.sensitivityElastic/temp1Max;
                temp2Max = max(max(abs(DV.sensitivityHeat)));
                DV.sensitivityHeat = DV.sensitivityHeat/temp2Max;
                
                g1 = config.w1*DV.sensitivityElastic+config.w2*DV.sensitivityHeat; % Calculate the weighted volume fraction change sensitivity.
            else
                g1 = config.w1*DV.sensitivityElastic;
            end
            
            % Filter the g1 sensitivies
            [g1]   = DV.check( config.nelx, config.nely,config.rmin,DV.x,g1);
            G1 = g1 - DV.lambda1 +1/(DV.mu1)*( targetFraction_v1-fractionCurrent_V1Local); % add in the lagrangian
            DV.w = DV.w+config.timestep*G1; % update the volume fraction.
            DV.w = max(min( DV.w,1),0);    % Don't allow the    vol fraction to go above 1 or below 0
            DV.lambda1 =  DV.lambda1 -1/(DV.mu1)*(targetFraction_v1-fractionCurrent_V1Local)*config.volFractionDamping;
        end
        
        % ----------------------------------
        % ORTHO DISTRIBUTION OPTIMIZATION
        % ----------------------------------
        %         function [] = OptimizeOrthoDistribution(obj,DV,config, matProp, masterloop)
        %             DV = DV.CalculateOthogonalDistributionSensitivity(config, matProp, masterloop);
        %             DV.sensitivityElastic = check( config.nelx, config.nely,config.rmin,DV.x,DV.sensitivityElastic);
        %             % move= 0.1* 20/(20+masterloop);
        %             move = config.orthDistMoveLimit;
        %             config.orthDistMoveLimit= config.orthDistMoveLimit* 10/(10+masterloop);
        %             %-----------------------
        %             %
        %             % Update design var.
        %             %-----------------------
        %             for ely = 1:config.nely
        %                 for elx = 1:config.nelx
        %                     if(DV.sensitivityElastic(ely,elx)<0.05)
        %                         DV.d(ely,elx) =  max(  DV.d(ely,elx)-move,config.minDorth);
        %                     end
        %
        %                     if(DV.sensitivityElastic(ely,elx)>0.05)
        %                         DV.d(ely,elx) =  min(  DV.d(ely,elx)+ move,config.maxDorth);
        %                     end
        %
        %                 end
        %             end
        %         end
        
        % ----------------------------------
        % ROTATION OPTIMIZATION
        % ----------------------------------
        function DV = OptimizeRotation(obj,DV,config, matProp, masterloop)
            %                 move= 0.1* 20/(20+masterloop);
            % allow multiple loading cases.
            [~, t2] = size(config.loadingCase);
            
            epsilon = pi/180; % 1 DEGREES ACCURACY
            elementsInRow = config.nelx+1;
            
            for ely = 1:config.nely
                rowMultiplier = ely-1;
                for elx = 1:config.nelx
                    rhoSIMP =  DV.x(ely,elx);
                    if(rhoSIMP>config.noNewMesoDesignDensityCutOff)
                        
                        % -------------------
                        % STEP 1, GET THE DISPLACEMENT FOR THIS NODE
                        % -------------------
                        nodes1=[rowMultiplier*elementsInRow+elx;
                            rowMultiplier*elementsInRow+elx+1;
                            (rowMultiplier +1)*elementsInRow+elx+1;
                            (rowMultiplier +1)*elementsInRow+elx];
                        
                        xNodes = nodes1*2-1;
                        yNodes = nodes1*2;
                        NodeNumbers = [xNodes(1) yNodes(1) xNodes(2) yNodes(2) xNodes(3) yNodes(3) xNodes(4) yNodes(4)];
                        UallCaseForElement = DV.U(1:t2,NodeNumbers);
                        U = UallCaseForElement;
                        
                        % -------------------
                        % STEP 2, SET UP GOLDEN RATIO METHOD TO FIND
                        % OPTIMAL THETA FOR ROTATION
                        % -------------------
                        
                        n = 0;
                        x0 = config.minRotation; %lower_bracket;
                        x3 = config.maxRotation;% higher_bracket;
                        leng = x3-x0;
                        grleng = leng*config.gr ; % golden ratio lenth
                        x1 = x3 - grleng;
                        x2 = x0 + grleng;
                        rhoSIMP =  DV.x(ely,elx);
                        mat1Frac  =[];% DV.w(ely,elx);
                        Exx = DV.Exx(ely,elx);
                        Eyy = DV.Eyy(ely,elx);
                        
                        thetaSubSystem = DV.thetaSub(ely,elx);
                        penaltyValue=DV.penaltyTheta(ely,elx);
                        lagraMultiplier=DV.lambdaTheta(ely,elx);
                        
                        %                         orthD = DV.d(ely,elx);
                        %fx1 = obj.EvaluteARotation(U,rhoSIMP, mat1Frac,Exx,Eyy,x1,thetaSubSystem,penaltyValue,lagraMultiplier,matProp, config,DV.maxElemStraniEnergy);
                        
                        
                        fx1= obj.EvaluteARotation(U,rhoSIMP, mat1Frac,Exx,Eyy,x1,thetaSubSystem,penaltyValue,lagraMultiplier,matProp, config,DV.maxElemStraniEnergy);
                        fx2 = obj.EvaluteARotation(U,rhoSIMP, mat1Frac,Exx,Eyy,x2,thetaSubSystem,penaltyValue,lagraMultiplier,matProp, config,DV.maxElemStraniEnergy);
                        
                        %                         if(masterloop>5)
                        %                             debug = 1;
                        %                         else
                        debug=0;
                        %                         end
                        verbosity = 0;
                        
                        if(   debug == 1)
                            xtemp = x0:pi/180:x3;
                            ytemp = zeros(1, size(xtemp,2));
                            count = 1;
                            for thetaTemp = xtemp
                                ytemp(count)= obj.EvaluteARotation(U,rhoSIMP, mat1Frac,Exx,Eyy,thetaTemp,thetaSubSystem,penaltyValue,lagraMultiplier,matProp, config,DV.maxElemStraniEnergy);
                                count = count+1;
                            end
                            figure(2)
                            subSysXvalus = [x0 DV.thetaSub(ely,elx) x3];
                            subSysYvalus = [min(ytemp) max(ytemp) max(ytemp)];
                            plot(xtemp,ytemp);
                            hold on
                            stairs(subSysXvalus,subSysYvalus)
                            hold off
                            title(sprintf('Lagrangian Function for Element x = %i, y = %i',elx,ely));
                            nothin = 1;
                        end
                        
                        
                        while(1 == 1)
                            if(debug == 1 && verbosity ==1)
                                str = sprintf('loop# = %d, x0 = %f, x1 = %f, x2 = %f, x3 = %f, fx1 = %f, fx2 = %f\n', n, x0, x1, x2, x3, fx1, fx2); display(str);
                            end
                            
                            if(fx1<=fx2) % less than or equal
                                % x0 = x0; % x0 stays the same
                                x3 = x2; % the old x2 is now x3
                                x2 = x1; % the old x1 is now x2
                                fx2 = fx1;
                                leng = x3 - x0; % find the length of the interval
                                x1 = x3 - leng*config.gr; % find golden ratio of length, subtract it from the x3 value
                                fx1 = obj.EvaluteARotation(U,rhoSIMP, mat1Frac,Exx,Eyy,x1,thetaSubSystem,penaltyValue,lagraMultiplier,matProp, config,DV.maxElemStraniEnergy);% calculate the fx
                                
                            elseif(fx1>fx2) % greater than
                                x0 = x1; % the old x1 is now x0
                                x1 = x2; % the old x2 is now the new x1
                                fx1 = fx2;
                                % x3 = x3; % x3 stays the same.
                                
                                leng = (x3 - x0); % find the length of the interval
                                x2 = x0 + leng*config.gr; % find golden ratio of length, subtract it from the x3 value
                                fx2 = obj.EvaluteARotation(U,rhoSIMP, mat1Frac,Exx,Eyy,x2,thetaSubSystem,penaltyValue,lagraMultiplier,matProp, config,DV.maxElemStraniEnergy);  % calculate the fx
                            end
                            
                            % check to see if we are as close as we want
                            if(leng < epsilon || n>100)
                                break;
                            end
                            n = n +1; % increment
                            
                        end
                        
                        % -------------------
                        % STEP 3, RECORD THE OPTIMAL THETA
                        % -------------------
                        minTvalue = (x2 + x3)/2;
                        moveLimit = config.rotationMoveLimit;
                        
                        % max move limit =  half the diff to optimal
                        diffT = abs(minTvalue-DV.t(ely,elx));
                        moveLimit=min(moveLimit,diffT*0.1);
                        
                        if(minTvalue>DV.t(ely,elx)+moveLimit)
                            DV.t(ely,elx)= DV.t(ely,elx)+moveLimit;
                        elseif(minTvalue<DV.t(ely,elx)-moveLimit)
                            DV.t(ely,elx)= DV.t(ely,elx)-moveLimit;
                        else
                            DV.t(ely,elx)=minTvalue;
                        end
                    end
                end
            end
        end
        
        % ---------------------------
        % EVALUTE THE OBJECTIVE FUNCTION FOR A ROTATION
        %----------------------------
        function lagrangianValue = EvaluteARotation(~,U,topDensity, material1Fraction,Exx,Eyy,thetaSys,thetaSubSystem,penaltyValue,lagraMultiplier,matProp, config,maxElemStraniEnergy)
            K = matProp.getKMatrixTopExxYyyRotVars(config,topDensity,Exx, Eyy,thetaSys,material1Fraction);
            % LOOP OVER LOADING CASES.
            % U'S ROWS ARE UNIQUE LOADING CASES
            % EACH ROW CONTAINS 8 VALUES FOR THE 8 DOF OF THE ELEMENT
            % allow multiple loading cases.
            [~, t2] = size(config.loadingCase);
            term1=0;
            for i = 1:t2
                Ucase = U(i,:)';
                term1= term1+Ucase'*K*Ucase;
            end
            term1=-term1;
            %             term1=-term1/maxElemStraniEnergy;
            
            %             term2 = penaltyValue/2*(thetaSys-thetaSubSystem)^2;
            term2 = penaltyValue*(thetaSys-thetaSubSystem)^2;
            %             term3 = lagraMultiplier*(thetaSys-thetaSubSystem);
            %             lagrangianValue=term1+term2+term3;
            %             normalizer=penaltyValue/2*(pi/4)^2;
            %                term2=term2/normalizer;
            
            %                lagrangianValue=term1+term2+term3;
            lagrangianValue=term1+term2;
            
            
        end
        
        
        % ----------------------------------
        % E_xx and E_yy  OPTIMIZATION
        % ----------------------------------
        function [DV] = OptimizeExxEyy(obj,DV,config, matProp, masterloop)
           DV= OptimizeExxEyy_V3(obj,DV,config, matProp, masterloop);
        end
        
        
        % Version 2
        function [DV] = OptimizeExxEyy_V2(obj,DV,config, matProp, masterloop)
            DV = DV.CalculateExxEyySensitivity(config, matProp, masterloop);
            DV.sensitivityElastic = DV.check( config.nelx, config.nely,config.rmin,DV.x,DV.sensitivityElastic);
            DV.sensitivityElasticPart2 = DV.check( config.nelx, config.nely,config.rmin,DV.x,DV.sensitivityElasticPart2);
            
            
            if(config.testingVerGradMaterail ==1)
                avgSensitivy = 0.5*( DV.sensitivityElastic+  DV.sensitivityElasticPart2);
                DV.sensitivityElastic =avgSensitivy;
                DV.sensitivityElasticPart2 =avgSensitivy;
            end
            
            
            %-----------------------
            %
            % Update design var.
            %-----------------------
            largest=1e8;
            divideByLargestFlag=0;
            l1 = 0; l2 = largest;% move = 0.2;
            
            E_target =(config.v1*matProp.E_material1+config.v2*matProp.E_material2)/(config.v1+config.v2);
            DV.targetAverageE = E_target;
            
            move = matProp.E_material1*0.05;
            minimum = matProp.E_material2;
            
            % ----------------
            % Exx
            % ----------------
            ExxNew = DV.Exx;
            EyyNew = DV.Eyy;
            
            
            offsetup = 10000;
            
            
            totalMaterial = sum(sum(DV.x));
            
            term1Exx = DV.sensitivityElastic;
            term1Eyy= DV.sensitivityElasticPart2;
            
            useOptimalCriteria = 0;
            
            % ---------------------------------------------------
            %
            % TARGET E AS THE CONSTRAINT
            %
            % ---------------------------------------------------
            if(config.useTargetMesoDensity~=1)
                while (l2-l1 > 1e-4)
                    lmid = 0.5*(l2+l1);
                    
                    
                    %---------------------------
                    % Add ATC terms
                    %
                    % Add in term 2 and 3 of the numerator for the consistency
                    % constraints for ATC optimization.
                    %---------------------------
                    term2Exx = DV.penaltyExx.*(ExxNew- DV.ExxSub);
                    term2Eyy = DV.penaltyEyy.*(EyyNew- DV.EyySub);
                    
                    term3Exx = DV.lambdaExx;
                    term3Eyy = DV.lambdaEyy;
                    
                    numeratorExx = term1Exx+term2Exx+term3Exx;
                    numeratorEyy = term1Eyy+term2Eyy+term3Eyy;
                    
                    % Don't allow negative
                    min1= min(min(numeratorExx));
                    min2= min(min(numeratorEyy));
                    min3 = min(min1,min2);
                    if(min1<=0 || min2<=0)
                        numeratorExx = numeratorExx-min3+1;
                        numeratorEyy = numeratorEyy-min3+1;
                    end
                    
                    % scale the sensitivies to make them easiler to work with if
                    % they are small.
                    max1= max(max(abs(numeratorExx)));
                    max2= max(max(abs(numeratorEyy)));
                    if(max1<=100 || max2<=100)
                        numeratorExx = numeratorExx*offsetup;
                        numeratorEyy = numeratorEyy*offsetup;
                    end
                    
                    if(max1>=10000000 || max2<=10000000)
                        numeratorExx = numeratorExx/offsetup;
                        numeratorEyy = numeratorEyy/offsetup;
                    end
                    
                    %                     ExxNew = max( minimum - EyyNew,  max(DV.Exx-move ,  min(  min(DV.Exx.*sqrt(DV.sensitivityElastic     ./lmid),DV.Exx+move ),matProp.E_material1)));
                    %                     EyyNew = max(minimum -  ExxNew,  max(DV.Eyy-move ,  min(  min(DV.Eyy.*sqrt(DV.sensitivityElasticPart2./lmid),DV.Eyy+move ),matProp.E_material1)));
                    ExxNew = max( minimum - EyyNew,  max(DV.Exx-move ,  min(  min(DV.Exx.*sqrt(numeratorExx     ./lmid),DV.Exx+move ),matProp.E_material1)));
                    EyyNew = max(minimum -  ExxNew,  max(DV.Eyy-move ,  min(  min(DV.Eyy.*sqrt(numeratorEyy./lmid),DV.Eyy+move ),matProp.E_material1)));
                    
                    
                    totalExx =DV.x.*ExxNew;
                    totalEyy = DV.x.* EyyNew;
                    avgE = (totalExx+totalEyy)/2;
                    averageElasticLocal= sum(sum(avgE))/totalMaterial;
                    %               averageElasticLocal = (sum(sum(EyyNew.*Xtemp))+sum(sum(ExxNew.*Xtemp)))/neSolid;
                    %               averageElasticLocal=averageElasticLocal/2; % Becuse Eyy and Exx are from one element, so to get the average divide by 2
                    if E_target- averageElasticLocal<0;
                        l1 = lmid;
                    else
                        l2 = lmid;
                    end
                end
            else
                
                %------------------------------
                % Use lagrangian update scheme
                %------------------------------
                ExxNew=DV.Exx;
                EyyNew=DV.Eyy;
                DV.penaltyMesoDensity=1;
                
                % just loop 20 times for now.
                lagrangianArray = [];
                avergeLagrangians=0;
                OldavergeLagrangians=0;
                
                minTerm1 = min(min(min(term1Exx)),min(min(term1Eyy)));
                maxTerm1 = max(max(max(term1Exx)),max(max(term1Eyy)));
                rangeTerm1 = maxTerm1-minTerm1;
                
                term1Exx=1/rangeTerm1*term1Exx;
                term1Eyy=term1Eyy/rangeTerm1;
                
                l1 = -1; l2 = 1;% move = 0.2;
                sumDensity =0;
                %    while (l2-l1 > 1e-4)
                
                v_total = config.nelx*config.nely*config.totalVolume;
                for i = 1:20000
                    %                     lmid = 0.5*(l2+l1);
                    
                    [dDensityEyy, dDensityExx,rhoValue] = obj.CalculateDensitySensitivityandRho(ExxNew/matProp.E_material1,EyyNew/matProp.E_material1,DV.t,DV.ResponseSurfaceCoefficents,config,matProp);
                    
                    rhoValue=max(0.05,min(rhoValue,1));
                    
                    dDensityExx= ones(size(dDensityExx));
                    dDensityEyy= ones(size(dDensityExx));
                    
                    mesoDensity=sum(sum(rhoValue.*DV.x))/(config.nelx*config.nely*config.totalVolume);
                    diffMesoDensity = mesoDensity-config.targetExxEyyDensity;
                    %                         diffMesoDensity = config.targetExxEyyDensity-mesoDensity;
                    
                    term2_b_Exx=(dDensityExx.*DV.x)/v_total;
                    term2_b_Eyy=(dDensityEyy.*DV.x)/v_total;
                    
                    term2A=DV.penaltyMesoDensity*diffMesoDensity;
                    
                    
                    term2Exx=term2A*term2_b_Exx;
                    term2Eyy=term2A*term2_b_Eyy;
                    
                    %
                    term3Exx=DV.lagrangianMesoDensity*term2_b_Exx;
                    term3Eyy=DV.lagrangianMesoDensity*term2_b_Eyy;
                    
                    term4Exx = DV.penaltyExx.*(ExxNew- DV.ExxSub)/matProp.E_material1;
                    term4Eyy = DV.penaltyEyy.*(EyyNew- DV.EyySub)/matProp.E_material1;
                    
                    term5Exx = DV.lambdaExx;
                    term5Eyy = DV.lambdaEyy;
                    
                    %                         term1Multiplier=1;
                    %
                    %                                             d_lagrangianXX=term1Exx +term2Exx+term3Exx+term4Exx+term5Exx;
                    %                                             d_lagrangianYY=term1Eyy +term2Eyy+term3Eyy+term4Eyy+term5Eyy;
                    
                    d_lagrangianXX=term1Exx +term2Exx+term3Exx+term4Exx;
                    d_lagrangianYY=term1Eyy +term2Eyy+term3Eyy+term4Eyy;
                    
                    %                     d_lagrangianXX=term1Exx +lmid+term4Exx+term5Exx;
                    %                     d_lagrangianYY=term1Eyy +lmid+term4Eyy+term5Eyy;
                    
                    %                         d_lagrangianXX=term3Exx;
                    %                         d_lagrangianYY=term3Eyy;
                    
                    deltaT=matProp.E_material1/10;
                    ddOffsetXX=deltaT*d_lagrangianXX;
                    ExxNew2=DV.Exx+ddOffsetXX;
                    
                    ddOffsetYY = deltaT*d_lagrangianYY;
                    EyyNew2=DV.Eyy+ddOffsetYY;
                    
                    ExxNew = max(0,max( minimum - EyyNew2,  max(DV.Exx-move ,  min(  min(ExxNew2,DV.Exx+move ),matProp.E_material1))));
                    EyyNew = max(0,max(minimum -  ExxNew2,  max(DV.Eyy-move ,  min(  min(EyyNew2,DV.Eyy+move ),matProp.E_material1))));
                    %                         ExxNew =max(0,min(ExxNew2,matProp.E_material1));
                    %                         EyyNew=max(0,min(EyyNew2,matProp.E_material1));
                    %                           DV.lagrangianMesoDensity = DV.lagrangianMesoDensity+term2A;
                    %                     multiplier=max(DV.lagrangianMesoDensity/100,1);
                    lagrangianArray=[lagrangianArray DV.lagrangianMesoDensity];
                    %                         DV.lagrangianMesoDensity = DV.lagrangianMesoDensity-sign(term2A)* term2A*multiplier;
                    %                     DV.lagrangianMesoDensity = DV.lagrangianMesoDensity-diffMesoDensity*multiplier;
                    DV.lagrangianMesoDensity = DV.lagrangianMesoDensity-diffMesoDensity*DV.penaltyMesoDensity;
                    
                end
                
                text1 =    sprintf('\nExxNew\t\t\t\t\t[1,1   5,1      5,5],  %f ,%f %f',ExxNew(1,1) ,ExxNew(5,1),ExxNew(5,5));
                text2 =    sprintf('ExxSub\t\t\t\t\t[1,1   5,1      5,5],  %f ,%f %f',DV.ExxSub(1,1),DV.ExxSub(5,1),DV.ExxSub(5,5));
                text3=    sprintf('Exx\t\t\t\t\t\t[1,1   5,1      5,5],  %f ,%f %f',DV.Exx(1,1),DV.Exx(5,1),DV.Exx(5,5));
                text33 =    sprintf('rhoValue \t\t\t\t[1,1   5,1      5,5],  %f ,%f %f',rhoValue(1,1),rhoValue(5,1),rhoValue(5,5));
                
                text4 =    sprintf('DiffXXNew \t\t\t\t[1,1   5,1      5,5],  %f ,%f %f',ExxNew(1,1) - DV.ExxSub(1,1),ExxNew(5,1) - DV.ExxSub(5,1),ExxNew(5,5) - DV.ExxSub(5,5));
                %                     text1 =    sprintf('numeratorExx\t\t\t[1,1   5,1      5,5],  %f ,%f, %f',numeratorExx(1,1),numeratorExx(5,1),numeratorExx(5,5));
                %                     text2 =    sprintf('combinedTermsExx/lmid\t[1,1   5,1      5,5],  %f ,%f %f\n',combinedTermsExx(1,1)/lmid,combinedTermsExx(5,1)/lmid,combinedTermsExx(5,5)/lmid);
                text5 =    sprintf('DiffXX \t\t\t\t\t[1,1   5,1      5,5],  %f ,%f %f',DV.Exx(1,1) - DV.ExxSub(1,1),DV.Exx(5,1) - DV.ExxSub(5,1),DV.Exx(5,5) - DV.ExxSub(5,5));
                %                     text6 =    sprintf('lmid\t\t\t\t\t\t\t%f',                      lmid);
                text7 =    sprintf('Meso Density [Target, Current, Diff] \t\t\t%f\t%f\t%f',  config.targetExxEyyDensity,mesoDensity,diffMesoDensity);
                %   text8 =    sprintf('Loop count\t\t\t%f', i );
                
                disp(text1)
                disp(text2)
                disp(text3)
                disp(text33)
                disp(text4)
                disp(text5)
                %                     disp(text6)
                disp(text7)
                %   disp(text8);
                
                % -----------------------
                % Set the valeus.
                % -----------------------
                DV.Exx =sqrt( ExxNew./  DV.Exx);
                
                DV.Eyy =sqrt( EyyNew./  DV.Eyy );
                
                %                       averageExx = mean(mean(ExxNew))
                %                      averageEyy= mean(mean(EyyNew))
                
                
            end
            
            if(config.testingVerGradMaterail ==1)
                averageNewE = 0.5*(ExxNew+EyyNew);
                ExxNew=averageNewE;
                EyyNew=averageNewE;
            end
            
            
            DV.Exx=ExxNew ;
            DV.Eyy=EyyNew ;
            
            
        end
        
        % ----------------------------------
        % E_xx and E_yy  OPTIMIZATION
        %
        %   Version 3
        % ----------------------------------
        function [DV] = OptimizeExxEyy_V3(obj,DV,config, matProp, masterloop)
            DV = DV.CalculateExxEyySensitivity(config, matProp, masterloop);
            DV.sensitivityElastic = DV.check( config.nelx, config.nely,config.rmin,DV.x,DV.sensitivityElastic);
            DV.sensitivityElasticPart2 = DV.check( config.nelx, config.nely,config.rmin,DV.x,DV.sensitivityElasticPart2);
            
          
             if(config.macro_meso_iteration>=2 && mod(masterloop,3)==1)
                 deltaT=2;
                   diffExx = DV.ExxSub-DV.Exx;
                 diffEyy = DV.EyySub-DV.Eyy  ;
                  
                 DV.lambdaExx= DV.lambdaExx+deltaT *diffExx;
                DV.lambdaEyy=  DV.lambdaEyy+deltaT*diffEyy;
                disp('Updated Lambda Values Exx Eyy')
                
            end
            
            
            %-----------------------
            %
            % Update design var.
            %-----------------------
            largest=1e8;
          
            l1 = 0; l2 = largest;% move = 0.2;
            
          
            
            move = matProp.E_material1*0.05;
            minimum = matProp.E_material2;
            
            % ----------------
            % Exx
            % ----------------
            ExxNew = DV.Exx;
            EyyNew = DV.Eyy;
            
            totalMaterial = sum(sum(DV.x));
            
            term1Exx = DV.sensitivityElastic;
            term1Eyy= DV.sensitivityElasticPart2;
            
            smallestLambdExx = min(min(DV.lambdaExx));
           smallestLambdEyy = min(min(DV.lambdaEyy));
           smallestOfTwo = min(smallestLambdExx,smallestLambdEyy);
            
            term2Exx =( DV.lambdaExx-smallestOfTwo).*DV.penaltyExx;
            term2Eyy = (DV.lambdaEyy-smallestOfTwo).*DV.penaltyEyy;
                     
            
            % ---------------------------------------------------
            %
            % TARGET AVG MESO DENSITY AS CONSTRAINT.
            % Update 3 (idea) 
            %
            % ---------------------------------------------------
            negValueFlag=0;
            negValueAdder=0;
            l1 = 0; l2 = largest;% move = 0.2;
            sumDensity =0;
            while (l2-l1 > 1e-4)
                lambda1 = 0.5*(l2+l1);
%                 [dDensityEyy, dDensityExx,rhoValue] = obj.CalculateDensitySensitivityandRho(ExxNew/matProp.E_material1,EyyNew/matProp.E_material1,DV.t,DV.ResponseSurfaceCoefficents,config,matProp);
                
                % testing. Set equal to one for now.
                dDensityExx=ones(size(term1Exx));
                dDensityEyy=ones(size(term1Exx));
                
               
                combinedTermsExx=(term1Exx+term2Exx)./(lambda1*dDensityExx);
                combinedTermsEyy=(term1Eyy+term2Eyy)./(lambda1*dDensityEyy);
                
                
%                 targetExx = ExxNew.*combinedTermsExx;
%                 targetEyy = EyyNew.*combinedTermsEyy;
                 targetExx = DV.Exx.*combinedTermsExx;
                targetEyy = DV.Eyy.*combinedTermsEyy;
                ExxNew = max(0,max( minimum - EyyNew,  max(DV.Exx-move ,  min(  min(targetExx,DV.Exx+move ),matProp.E_material1))));
                EyyNew = max(0,max(minimum -  ExxNew,  max(DV.Eyy-move ,  min(  min(targetEyy,DV.Eyy+move ),matProp.E_material1))));
                
                
                
                for i = 1:config.nelx
                    for j = 1:config.nely
                        % scale down the X and Y
                        x=ExxNew(j,i)/matProp.E_material1;
                        y=EyyNew(j,i)/matProp.E_material1;
                        theta=DV.t(j,i);
                        [~, ~,estimateElementDensity] = obj.CalculateDensitySensitivityandRho(x,y,theta,DV.ResponseSurfaceCoefficents,config,matProp);
                        
                        estimateElementDensity= min(max(estimateElementDensity,0.05),1);%1 is max, 0.5 is min
                        eleDensity = DV.x(j,i)*estimateElementDensity;
                        sumDensity =sumDensity+eleDensity;
                        
                        
                    end
                end
                sumDensity = sumDensity/(config.nelx*config.nely*config.totalVolume);
                
                if config.targetExxEyyDensity- sumDensity<0;
                    l1 = lambda1;
                    %                          l2 = lmid;
                else
                    l2 = lambda1;
                    %                         l1 = lmid;                   
                    
                end
            end
            multiplier= 10000;
            text1 =    sprintf('\nExxNew\t\t\t\t\t[1,1   5,1      5,5],  %f ,%f %f',ExxNew(1,1) ,ExxNew(5,1),ExxNew(5,5));
            text2 =    sprintf('ExxSub\t\t\t\t\t[1,1   5,1      5,5],  %f ,%f %f',DV.ExxSub(1,1),DV.ExxSub(5,1),DV.ExxSub(5,5));
            text3=    sprintf('Exx\t\t\t\t\t\t[1,1   5,1      5,5],  %f ,%f %f',DV.Exx(1,1),DV.Exx(5,1),DV.Exx(5,5));
            text4 =    sprintf('DiffXXNew \t\t\t\t[1,1   5,1      5,5],  %f ,%f %f',ExxNew(1,1) - DV.ExxSub(1,1),ExxNew(5,1) - DV.ExxSub(5,1),ExxNew(5,5) - DV.ExxSub(5,5));
            text5 =    sprintf('combinedTermsExx\t[1,1   5,1      5,5],  %f ,%f %f',combinedTermsExx(1,1),combinedTermsExx(5,1),combinedTermsExx(5,5));
            text6 =    sprintf('penaltyExx*10000 \t\t\t\t\t[1,1   5,1      5,5],  %f ,%f %f',DV.penaltyExx(1,1)*multiplier,DV.penaltyExx(5,1)*multiplier,DV.penaltyExx(5,5)*multiplier);
            
            text7 =    sprintf('lambdaExx \t\t\t\t\t[1,1   5,1      5,5],  %f ,%f %f',DV.lambdaExx(1,1),DV.lambdaExx(5,1),DV.lambdaExx(5,5));
            
            text8 =    sprintf('lambda1 and density\t\t\t[  %f ,%f\n',lambda1,sumDensity);
            
            disp(text1)
            disp(text2)
            disp(text3)
            disp(text4)
            disp(text5)
            disp(text6)
            disp(text7)
            disp(text8)
            
            
            
            debug = 0;
            if(debug ==1)
                figure(2)
                p = plotResults;
                xplots=3;
                yplots =3;
                plotNum=1;
                
                subplot(xplots,yplots,plotNum);
                p. PlotArrayGeneric(dDensityEyy,'dDensityEyy');
                plotNum=plotNum+1;
                
                subplot(xplots,yplots,plotNum);
                p. PlotArrayGeneric(dDensityExx,'dDensityExx');
                plotNum=plotNum+1;
                
                subplot(xplots,yplots,plotNum);
                p. PlotArrayGeneric(numeratorExx,'completeExx');
                plotNum=plotNum+1;
                
                subplot(xplots,yplots,plotNum);
                p. PlotArrayGeneric(numeratorEyy,'completeEyy');
                plotNum=plotNum+1;
                
                subplot(xplots,yplots,plotNum);
                p. PlotArrayGeneric(combinedTermsExx/lmid,'combinedTermsExx/lmid');
                plotNum=plotNum+1;
                
                subplot(xplots,yplots,plotNum);
                p. PlotArrayGeneric(combinedTermsEyy/lmid,'combinedTermsEyy/lmid');
                plotNum=plotNum+1;
                
                subplot(xplots,yplots,plotNum);
                p. PlotArrayGeneric(DV.t-DV.thetaSub,'theta diff');
                plotNum=plotNum+1;
                
                subplot(xplots,yplots,plotNum);
                p. PlotArrayGeneric(DV.Eyy - DV.EyySub,'Eyy diff');
                plotNum=plotNum+1;
                
                subplot(xplots,yplots,plotNum);
                p. PlotArrayGeneric(DV.Exx-DV.ExxSub,'Exx diff');
                plotNum=plotNum+1;
                
                
            end
            
                          
                       
            % -----------------------
            % Set the valeus.
            % -----------------------
            DV.Exx =DV.Exx.*sqrt( ExxNew./  DV.Exx);
            
            DV.Eyy = DV.Eyy.*sqrt( EyyNew./  DV.Eyy );
          
            
            
        end
        
        
        function [EyySensitivty, ExxSensitivity,rhoValue] = CalculateDensitySensitivityandRho(obj,Exx,Eyy,theta,Coefficents,config,matProp)
            co = Coefficents;
            if(config.useThetaInSurfaceFit==1)
                
                % make it so that Exx is always larger
                temp1=Eyy;
                valueConditionTrue = Eyy>Exx;
                Eyy(valueConditionTrue)=Exx(valueConditionTrue);
                Exx(valueConditionTrue)=temp1(valueConditionTrue);
                %                 if(Eyy>Exx)
                %                     Exx=Eyy;
                %                     Eyy=Exx;
                %                 end
                
                
                
                % rhoValue= x(1)  + x(2)* exp(E_xx)  + x(3)* exp(E_yy)+x(4) *exp(theta) +x(5)*E_xx  + x(6)* E_yy +x(7)*theta+ x(8)*E_xx.*E_yy;
                % ExxSensitivity=  x(2)* exp(E_xx)  +x(5) + x(8)*E_yy;
                % EyySensitivty= x(3)* exp(E_yy) + x(6)+ x(8)*E_xx;
                rhoValue= co(1)+co(2)*Exx+co(3)*Eyy+co(4)*theta+co(5)*Exx.^2+co(6)* Eyy.^2+co(7)*theta.^2+co(8)*Exx.*Eyy+co(9)*Eyy.*theta+co(10)*Exx.*theta;
                ExxSensitivity =co(2)+2*co(5)*Exx+co(8)*Eyy+co(10)*theta;
                EyySensitivty = co(3)+2*co(6)* Eyy+co(8)*Exx+co(9)*theta;
                
                % Scale Up
                %                 rhoValue=rhoValue*scaleUp;
                %                 ExxSensitivity=ExxSensitivity*scaleUp;
                %                 EyySensitivty=EyySensitivty*scaleUp;
                
                %                 rhoValue(rhoValue>1)=1;
                %                 rhoValue(rhoValue<0)=1;
                
                
                
            else
                % obj. ResponseSurfaceCoefficents=[ 1.0000000000463e-05 9.99988184437107e-06 9.9998491550433e-06 -3.40115537230351e-11 -5.52110060132392e-12 -3.81038581303971e-11];
                
                minAllowed = 0.01;
                EyySensitivty=max(co(3)+co(5).*Exx+2*co(6).*Eyy,minAllowed);
                ExxSensitivity=max(co(2)+ 2*co(4).*Exx+co(5).*Eyy,minAllowed);
                rhoValue=   co(1) + co(2)*Exx + co(3)*Eyy + co(4)*Exx.^2 + co(5)*Exx.*Eyy + co(6)*Eyy.^2;
            end
            
        end
        
        %-----------------------------------
        % Meso Optimization
        %-----------------------------------
        function [DVmeso] = MesoDensityOptimization(~,mesoConfig, DVmeso,old_muMatrix,penaltyValue,macroElemProps)
            ne = mesoConfig.nelx*mesoConfig.nely; % number of elements
            %               dH_total=[DVmeso.d11;
            %                     DVmeso.d12;
            %                     DVmeso.d22;
            %                     DVmeso.d33];
            Diff_Sys_Sub =  (macroElemProps.D_subSys- macroElemProps.D_sys);
            localD = zeros(3,3);
            for e = 1:ne
                
                [x,y]= DVmeso.GivenNodeNumberGetXY(e);
                xx=DVmeso.x(y,x); % =min(optimalEta, designVars.x+move)
                %                  term1 = 10*xx^9;
                %                  power = 1/4;
                %                  term1 = power*xx^(power-1);
                term1=2*xx;
                
                
                
                rowIndex = [1,1,2,3];
                columnIndex = [1,2,2,3];
                
                dH = zeros(3,3);
                dH(1,1) = DVmeso.d11(y,x);
                dH(1,2) = DVmeso.d12(y,x);
                dH(2,2) = DVmeso.d22(y,x);
                dH(3,3) = DVmeso.d33(y,x);
                
                localD(1,1) = DVmeso.De11(y,x);
                localD(1,2) = DVmeso.De11(y,x);
                localD(2,2) = DVmeso.De11(y,x);
                localD(3,3) = DVmeso.De11(y,x);
                
                Diff_Sys_Sub =  (localD- macroElemProps.D_sys);
                
                constraintCount = 0;
                term2=0;
                %                 term1=0;
                for k = [1 2 3 ]
                    %                     term1=  dH(1,1)+  dH(1,2)+  dH(2,2)+  dH(3,3);
                    i = rowIndex(k);
                    j = columnIndex(k);
                    Ctemp = dH(i,j)*(-old_muMatrix(i,j)-penaltyValue*Diff_Sys_Sub(i,j));
                    term2 =term2 +Ctemp;
                    constraintCount=constraintCount+1;
                end
                
                dL = term1+term2;
                delta = 0.1;
                optimalEta=xx+delta*dL;
                move = 0.02;
                DVmeso.x(y,x)=  max(0.01,max(xx-move,min(1.,min(xx+move,optimalEta))));
                
                DVmeso.x([10:13],[10:13])=1;
            end
        end
    end
end

