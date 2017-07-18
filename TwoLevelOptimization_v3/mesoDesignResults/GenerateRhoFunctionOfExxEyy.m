function []=GenerateRhoFunctionOfExxEyy(config)
close all
% p = plotResults;
DV = DesignVars(config);
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



% Get the Exx field
outname = sprintf('./out%i/ExxValues%i.csv',folderNum,macro_meso_iteration);
ExxMacro = csvread(outname);

% Get the Eyy field
outname = sprintf('./out%i/EyyValues%i.csv',folderNum,macro_meso_iteration);
EyyMacro =csvread(outname);

% Get the Theta field
outname = sprintf('./out%i/ThetaValues%i.csv',folderNum,macro_meso_iteration);
ThetaMacro = csvread(outname);



matProp=MaterialProperties;
DV.x = xxx;
rhoArray = [];
ExxArray = [];
EyyArray = [];
thetaArray=[];
MacroExxColumn=[];
MacroEyyColumn=[];
MacroThetaColumn=[];
RhoColumn=[];

if(config.validationModeOn==1)
    ne= config. validationGridSizeNelx ^3;
end
ne

testing=0;

if 1==1
    for e = 1:ne %ne:-1:1
        strangeResultsFlag=0;
        fprintf('element %i of %i\n',e,ne);
        macroElementProps.elementNumber=e;
        elementNumber=e;
        if(config.validationModeOn==0)
            % ---------------
            %    Multiscale topology optimization case.
            % ---------------
            results = elementXYposition(macroElementProps.elementNumber,:);
            macroElementProps.yPos = results(1);
            macroElementProps.xPos = results(2);
            
            
            macroElementProps.densitySIMP = xxx(macroElementProps.yPos,macroElementProps.xPos );
            ActualThetaValue = ThetaMacro(macroElementProps.yPos,macroElementProps.xPos );
            ActualExx = ExxMacro(macroElementProps.yPos,macroElementProps.xPos );
            ActualEyy = EyyMacro(macroElementProps.yPos,macroElementProps.xPos );
        elseif(config.validationModeOn==1)
            % ---------------
            %    Meso Validation Case
            % ---------------
            macroElementProps.yPos = e;
            macroElementProps.xPos = 1;
            macroElementProps.densitySIMP =1;
            ActualThetaValue = ThetaMacro(e );
            ActualExx = ExxMacro(e );
            ActualEyy = EyyMacro(e);
            
            numValues = config.validationGridSizeNelx ^3;
            %             DV.Exx = ones(1,numValues);
            %             DV.Eyy= ones(1,numValues);
            %             DV.t= ones(1,numValues);
            %             DV.w= ones(1,numValues);
            %             DV.x =ones(1,numValues);
            DV.Exx = ones(numValues,1);
            DV.Eyy=ones(numValues,1);
            DV.t= ones(numValues,1);
            DV.w= ones(numValues,1);
            DV.x =ones(numValues,1);
            
        end
        
        
        
        % ---------------------------------
        % If SIMP density is above minimum, then find the equivalent macro
        % properties from the D_meso matrix.
        % ---------------------------------
        if(macroElementProps.densitySIMP>config.voidMaterialDensityCutOff)
            % save the psuedo strain values
            %         outname = sprintf('./out%i/psuedostrain_Ite%i_forElement%i.csv',folderNum,macro_meso_iteration,elementNumber);
            %     p =  macroElementProperties.psuedoStrain;
            %               p =  csvread(outname);
            
            if(testing==0)
                % get the final volume
                outname = sprintf('./out%i/volumeUsed_Ite%i_forElement%i.csv',folderNum,macro_meso_iteration,elementNumber);
                %              v =    configMeso.totalVolume;
                %     csvwrite(outname,v);
                if exist(outname, 'file') ~= 2
                    fprintf('File does not exist. Retry\n');
                    combinedTopologyOptimization('1', '1', '1','100', int2str(elementNumber));

                end
                v =  csvread(outname);
                if(v>1)
                    message = 'Volume is greater than 1???'
                    break
                end

                if(v<0)
                    message = 'Volume is less than 0???'
                    break
                end

                if(isreal(v)~=1)
                    message = 'Volume is imaginary???'
                    break
                end
            else
                v=1;
            end
            
            if(testing==0)
                outname = sprintf('./out%i/Dmatrix_%i_forElement_%i.csv',folderNum,macro_meso_iteration,elementNumber);
            else
                outname = sprintf('./out%i/DsystemIter%i_Element_%i.csv',folderNum,macro_meso_iteration,elementNumber);
            end
            
            if exist(outname, 'file') ~= 2
                continue;
            end
            %          outname = sprintf('./out%i/DsystemIter%i_Element_%i.csv',folderNum,macro_meso_iteration,elementNumber);
            try
                Din = csvread(outname);
            catch
                str = sprintf('error reading a file\n'); display(str);
                continue
                
            end
            
            Dcalculated= matProp.getDmatMatrixTopExxYyyRotVars(config, macroElementProps.densitySIMP,ActualExx, ActualEyy,ActualThetaValue,1);
            
            % -------------------
            % STEP 2, SET UP GOLDEN RATIO METHOD TO FIND
            % OPTIMAL THETA FOR ROTATION
            % -------------------
            
            n = 0;
            epsilon = pi/180; % 1 DEGREES ACCURACY
            x0 = 0; %lower_bracket;
            x3 =pi/2;% higher_bracket;
            leng = x3-x0;
            grleng = leng*config.gr ; % golden ratio lenth
            x1 = x3 - grleng;
            x2 = x0 + grleng;
            Dout = matProp.rotateDmatrix(config,x1, Din);
            sumZeroTerms =abs( Dout(1,3))+abs(Dout(2,3))+abs(Dout(3,1))+abs(Dout(3,2));
            fx1=sumZeroTerms;
            
            Dout = matProp.rotateDmatrix(config,x2, Din);
            sumZeroTerms =abs( Dout(1,3))+abs(Dout(2,3))+abs(Dout(3,1))+abs(Dout(3,2));
            fx2=sumZeroTerms;
            
            verbosity = 0;
            recordFx=[];
            
            while(1 == 1)
                if(verbosity ==1)
                    str = sprintf('loop# = %d, x0 = %f, x1 = %f, x2 = %f, x3 = %f, fx1 = %f, fx2 = %f\n', n, x0, x1, x2, x3, fx1, fx2); display(str);
                end
                
                if(fx1<=fx2) % less than or equal
                    % x0 = x0; % x0 stays the same
                    x3 = x2; % the old x2 is now x3
                    x2 = x1; % the old x1 is now x2
                    fx2 = fx1;
                    leng = x3 - x0; % find the length of the interval
                    x1 = x3 - leng*config.gr; % find golden ratio of length, subtract it from the x3 value
                    
                    Dout = matProp.rotateDmatrix(config,x1, Din);
                    fx1 =abs( Dout(1,3))+abs(Dout(2,3))+abs(Dout(3,1))+abs(Dout(3,2));
                    % fx1 = obj.EvaluteARotation(U,topDensity, material1Fraction,Exx,Eyy,x1,matProp, config); % calculate the fx
                    
                elseif(fx1>fx2) % greater than
                    x0 = x1; % the old x1 is now x0
                    x1 = x2; % the old x2 is now the new x1
                    fx1 = fx2;
                    % x3 = x3; % x3 stays the same.
                    
                    leng = (x3 - x0); % find the length of the interval
                    x2 = x0 + leng*config.gr; % find golden ratio of length, subtract it from the x3 value
                    %                 fx2 = obj.EvaluteARotation(U,topDensity, material1Fraction,Exx,Eyy,x2,matProp, config);  % calculate the fx
                    Dout = matProp.rotateDmatrix(config,x2, Din);
                    fx2 =abs( Dout(1,3))+abs(Dout(2,3))+abs(Dout(3,1))+abs(Dout(3,2));
                end
                
                % check to see if we are as close as we want
                if(leng < epsilon || n>100)
                    break;
                end
                n = n +1; % increment
                
            end
            Theta = (x2 + x3)/2;
            
            
            % plot the domain and range of the Kappa function
            if(1==0)
                thetaValuesToTest=0:epsilon:pi/2;
                for i = thetaValuesToTest
                    Dout = matProp.rotateDmatrix(config,i, Din);
                    fxRecord =abs( Dout(1,3))+abs(Dout(2,3))+abs(Dout(3,1))+abs(Dout(3,2));
                    recordFx=[recordFx fxRecord];
                end
                
                xValues = [0 Theta pi/2];
                yValues = [ 0 1 1]*max(recordFx);
                plot(thetaValuesToTest,recordFx);
                hold on
                stairs(xValues,yValues);
                titleText = sprintf('Kappa function for element %i',e);
                title(titleText);
                hold off
            end
            
            
            Din=Din*1/(macroElementProps.densitySIMP^(config.penal));
            Dout = matProp.rotateDmatrix(config,Theta, Din);
            
            % Scale UP based on the rho (x) density
            % Scale down on the denominator
            denominator = 1-matProp.v^2;
            %          denominator = 1;
            Etemp1=Dout(1,1)*denominator;
            Etemp2=Dout(2,2)*denominator;
            Theta = pi/2-Theta;
            
            
            if(ActualExx>ActualEyy)
                Exx=max(Etemp1,Etemp2);
                Eyy=min(Etemp1,Etemp2);
            else %ActualExx<ActualEyy
                Exx=min(Etemp1,Etemp2);
                Eyy=max(Etemp1,Etemp2);
            end
            diffDs=Din-Dcalculated;
            
            if(config.validationModeOn==0)
                if(Exx>matProp.E_material1)
                    Exx=matProp.E_material1 ;
                    strangeResultsFlag=1;
                end
                
                if(Eyy>matProp.E_material1)
                    Eyy=matProp.E_material1 ;
                    strangeResultsFlag=1;
                end
            end
            
            
            
            
            % Also, there is some difficulty when actualTheta is 0 or pi/2
              diffTheta = abs( ActualThetaValue-Theta);
            
            if(diffTheta>(pi/2-epsilon))
                if(Theta>pi/4)
                    Theta=pi/2-Theta;
                else
                     Theta=pi/2-Theta;
%                       Theta=Theta-pi/2;
                end
                
                %             str = sprintf('Theta on wrong boundary. Switching values. ')
%                 diffTheta = abs( ActualThetaValue-Theta); % The new Diff Theta
            end
            
            
            % in the case  where, Exx = Eyy, then the material is basically
            % isotropic and rotating to find the orthotropic orientaiton will
            % not work. In this case, set the theta to the actualTheta
            % the criteria is that Exx and Eyy must be within 0.2% of each
            % other's value.
            if(abs(100*(Exx-Eyy)/Exx)<5)
                %            if(abs(100*(ActualExx-ActualEyy)/ActualEyy)<5)
                Theta=ActualThetaValue;
                %             str = sprintf('Exx = Eyy, setting theta to sys theta')
            end
            
            Exx=max(0,Exx);
            Eyy=max(0,Eyy);
            
            diffX = ActualExx-Exx;
            diffY = ActualEyy-Eyy;
             diffTheta =  ActualThetaValue-Theta;
            
            
%             relativeErrorDiffExx = abs(diffX)/max(ActualExx,5);
%             relativeErrorDiffyy = abs(diffY)/max(ActualEyy,5);
%             relativeErrorDiffTheta = abs(diffTheta)/max(ActualThetaValue,0.05);
            
             relativeErrorDiffExx = abs(diffX)/mean([ActualExx Exx]);
            relativeErrorDiffyy = abs(diffY)/mean([ActualEyy Eyy]);
            relativeErrorDiffTheta = abs(diffTheta)/mean([ActualThetaValue Theta]);
            
            maxError = 1;
            LargeErroFlag = 0;
            
            if(1==1)
                if(ActualExx<20000 && ActualEyy< 500 && ActualThetaValue>pi/4  ) % target the low densities
%                     if(ActualThetaValue>pi/4  && ActualEyy*0.25>ActualExx )
                        if(relativeErrorDiffExx>maxError)
                            fprintf('%i Exx Large Error: Target  = %f mesovalue = %f\n', e,ActualExx,Exx)
                            LargeErroFlag = 1;
                        end
                        if(relativeErrorDiffyy>maxError)
                            fprintf('%i Eyy Large Error: Target  = %f mesovalue = %f\n', e,ActualEyy,Eyy)
                            LargeErroFlag = 1;
                        end
                        if( relativeErrorDiffTheta>maxError)
                            fprintf('%i Theta Large Error: Target  = %f mesovalue = %f\n', e,ActualThetaValue,Theta)
                            LargeErroFlag = 1;
                        end

                        if( LargeErroFlag ==1)
                            fprintf('More Data: Target: Value: Relative Error\n')
                            fprintf('Exx %f %f %f\n',ActualExx,Exx,relativeErrorDiffExx)
                            fprintf('Eyy %f %f %f\n',ActualEyy,Eyy,relativeErrorDiffyy)
                            fprintf('Theta %f %f %f\n',ActualThetaValue,Theta,relativeErrorDiffTheta)
                            fprintf('rho = %f\n\n',v);
                            %: Targets %f %f %f, Meso %f %f %f, Rho=%f\n',ActualExx,ActualEyy,ActualThetaValue,Exx,Eyy,Theta,v)
                             plotAnElement( config.macro_meso_iteration,e)
                        end
                    
%                     end
                end
            end
            
            
            
            
            
            %DV.x already saved.
            DV.Exx(macroElementProps.yPos,macroElementProps.xPos)=Exx;
            DV.Eyy(macroElementProps.yPos,macroElementProps.xPos)=Eyy;
            DV.t(macroElementProps.yPos,macroElementProps.xPos)=Theta;
            DV.w(macroElementProps.yPos,macroElementProps.xPos)=v;
            
            
            %                 objectiveValue = ObjectiveCalculateEffectiveVars(x,DmatrixIN, matProp,config);
            
            if(strangeResultsFlag==0)
                % The Exx and Eyy values here are already correct to take
                % into account the SIMP, rho value. So, you do not need to
                % scale the ActualExx or ActualEyy to compare the two. 
                
                rhoArray = [rhoArray;v];
                ExxArray=[ExxArray;Exx];
                EyyArray=[EyyArray;Eyy];
                thetaArray=[thetaArray;Theta];
                
%                 MacroExxColumn=[MacroExxColumn;ActualExx*macroElementProps.densitySIMP^(config.penal)];
%                 MacroEyyColumn=[MacroEyyColumn;ActualEyy*macroElementProps.densitySIMP^(config.penal)];
                
                  MacroExxColumn=[MacroExxColumn;ActualExx];
                MacroEyyColumn=[MacroEyyColumn;ActualEyy];
                MacroThetaColumn=[MacroThetaColumn;ActualThetaValue];
                RhoColumn=[RhoColumn;v];
            end
            
            
            
        else
            DV.Exx(macroElementProps.yPos,macroElementProps.xPos)=ActualExx;
            DV.Eyy(macroElementProps.yPos,macroElementProps.xPos)=ActualEyy;
            DV.t(macroElementProps.yPos,macroElementProps.xPos)=ActualThetaValue;
            DV.w(macroElementProps.yPos,macroElementProps.xPos)=0;
            
        end
    end
    
      if(config.validationModeOn==1)
        PlotType = 'MesoValidation' ;
    else
        PlotType=sprintf('Results From Iter %i',config.macro_meso_iteration);
    end
    
    
    
    % --------------------
    %
    %   PLOT a comparison
    %
    % --------------------
    p = plotResults;
    size(ExxMacro)
    size(DV.Exx)
    diffExx = ExxMacro- DV.Exx;
    diffEyy = EyyMacro- DV.Eyy;
    diffTheta = ThetaMacro- DV.t;
    
    logic = DV.x>config.voidMaterialDensityCutOff;
    zeroArray = zeros(size(diffExx));
    diffExx(~logic)=zeroArray(~logic);
    diffEyy(~logic)=zeroArray(~logic);
    diffTheta(~logic)=zeroArray(~logic);
    
    
    
    xplots = 2;
    yplots = 2;
    c= 1;
    figure(1)
    subplot(xplots,yplots,c);c=c+1;
    
    p.PlotArrayGeneric( diffExx, 'Macro - Sub , diffExx')
    subplot(xplots,yplots,c);c=c+1;
    p.PlotArrayGeneric( diffEyy, 'Macro - Sub ,diffEyy')
    subplot(xplots,yplots,c);c=c+1;
    p.PlotArrayGeneric( diffTheta, 'Macro - Sub ,diffTheta')
    nameGraph = sprintf('./%s_MesoDesignExxEyyThetaActualVsTarget%i.png',PlotType, config.macro_meso_iteration);
    print(nameGraph,'-dpng');
    close all
    
    % ------------------
    % Plot the error as positions on the design
    % --------------
    diffExx=abs(diffExx);
    sumABSVAlues = abs(ExxMacro)+abs( DV.Exx);
    diffExx=diffExx./(sumABSVAlues./2); %
    
    diffEyy=abs(diffEyy);
    sumABSVAlues = abs(EyyMacro)+abs(DV.Eyy);
    diffEyy=diffEyy./(sumABSVAlues./2); %
    
    diffTheta=abs(diffTheta);
    sumABSVAlues = abs(ThetaMacro)+abs( DV.t);
    diffTheta=diffTheta./(sumABSVAlues./2); %
    
    totalError = diffExx+(diffEyy)+(diffTheta);
    
    
    
    xplots = 2;
    yplots = 2;
    c= 1;
    figure(1)
    subplot(xplots,yplots,c);c=c+1;
    
    p.PlotArrayGeneric( diffExx, 'Error Exx')
    subplot(xplots,yplots,c);c=c+1;
    p.PlotArrayGeneric( diffEyy, 'Error Eyy')
    subplot(xplots,yplots,c);c=c+1;
    p.PlotArrayGeneric( diffTheta, 'Error theta')
    subplot(xplots,yplots,c);c=c+1;
    p.PlotArrayGeneric( totalError, 'Error totalError')
    nameGraph = sprintf('./%s_MesoDesignExxEyyThetaActualVsTarget%i.png',PlotType, config.macro_meso_iteration);
    print(nameGraph,'-dpng');
    close all
    
%      Errorxx =abs( diffExx)./0.5*( ExxMacro+ DV.Exx);
%     Erroryy =abs( diffEyy)./0.5*( EyyMacro+ DV.Eyy);% EyyMacro- DV.Eyy;
%     ErrorTheta =abs( diffTheta)./0.5*( ThetaMacro+ DV.t);%  ThetaMacro- DV.t;
%     
%      xplots = 2;
%     yplots = 2;
%     c= 1;
%     figure(1)
%     subplot(xplots,yplots,c);c=c+1;
%     
%     p.PlotArrayGeneric( Errorxx, 'Error (diff)/avg, Exx')
%     subplot(xplots,yplots,c);c=c+1;
%     p.PlotArrayGeneric( Erroryy, 'Error (diff)/avg, Eyy')
%     subplot(xplots,yplots,c);c=c+1;
%     p.PlotArrayGeneric( ErrorTheta,  'Error (diff)/avg, Theta')
%     nameGraph = sprintf('./%s_MesoDesignExxEyyThetaErrorINPositions%i.png',PlotType, config.macro_meso_iteration);
%     print(nameGraph,'-dpng');
%     close all
%     
    
    % --------------------
    %
    %   SCALE the Sub Values
    %
    % --------------------
    
    if(config.ScaleTheSubSystemValuesToMeetVolumeConstraint==1)
    
        l1 = 0; l2 = 3;% move = 0.2;
        %             sumDensity =0;
        o=Optimizer;
        if(config.useTargetMesoDensity==1)
            target=config.targetExxEyyDensity;
            theta=DV.t;
        else
            target=config.targetAvgExxEyy;
            totalMaterial= sum(sum(DV.x));
        end
    
        fprintf('try scaling the Sub values\n');
    
        while (l2-l1 > 1e-6)
            lambda1 = 0.5*(l2+l1);
            ExxNew=DV.Exx*lambda1;
            EyyNew=DV.Eyy*lambda1;
    
            if(config.useTargetMesoDensity==1)
                [~, ~,rhoValue] = o.CalculateDensitySensitivityandRho(ExxNew/matProp.E_material1,EyyNew/matProp.E_material1,theta,DV.x,DV.ResponseSurfaceCoefficents,config,matProp,0);
                rhoValue=max(0,min(rhoValue,1));
                temp2 = sum(sum(rhoValue));
                sumDensity=temp2/(config.nelx*config.nely*config.totalVolume);
                currentValue=sumDensity;
            else
    
    
                totalExx =DV.x.*ExxNew;
                totalEyy = DV.x.* EyyNew;
                avgE = (totalExx+totalEyy)/2;
                averageElasticLocal= sum(sum(avgE))/totalMaterial;
    
                currentValue=averageElasticLocal;
            end
    
    
            fprintf('Target %f and current %f\n',target,currentValue);
            if target- currentValue<0;
                l2 = lambda1;
            else
                l1 = lambda1;
            end
        end
    
        DV.Exx=    DV.Exx*lambda1;
        DV.Eyy=      DV.Eyy*lambda1;
    
        fprintf('Final Lambda = %f with final value of %f\n',lambda1,currentValue);
    end
    
    % --------------------
    %
    %   SAVE the Sub Values and everything else
    %
    % --------------------
    
    %     folderCells = {sprintf('out%i',folderNum),'data'};
    if(config.validationModeOn==1)
        folderName='data';
    else
        folderName=sprintf('out%i',folderNum);
    end
    %     for i = 1:2
    %         folderName = char(folderCells(i));
    % save the Exx field
    %     outname = sprintf('./out%i/ExxSubSysValues%i.csv',folderNum,macro_meso_iteration);
    outname = sprintf('./%s/ExxSubSysValues%i.csv',folderName,macro_meso_iteration);
    csvwrite( outname,DV.Exx);
    
    % save the Eyy field
    %     outname = sprintf('./out%i/EyySubSysValues%i.csv',folderNum,macro_meso_iteration);
    outname = sprintf('./%s/EyySubSysValues%i.csv',folderName,macro_meso_iteration);
    csvwrite( outname,DV.Eyy);
    
    % save the Theta field
    %     outname = sprintf('./out%i/ThetaSubSysValues%i.csv',folderNum,macro_meso_iteration);
    outname = sprintf('./%s/ThetaSubSysValues%i.csv',folderName,macro_meso_iteration);
    csvwrite( outname, DV.t);
    
    % save the density field
    %     outname = sprintf('./out%i/densityUsedSubSysValues%i.csv',folderNum,macro_meso_iteration);
    outname = sprintf('./%s/densityUsedSubSysValues%i.csv',folderName,macro_meso_iteration);
    csvwrite( outname,  DV.w);
    
    %------------------------
    % Save the macro columns as well. This will help with future analysis
    % ----------------------------
    % save the MacroExxColumn
    %     outname = sprintf('./out%i/MacroExxColumn%i.csv',folderNum,macro_meso_iteration);
    outname = sprintf('./%s/MacroExxColumn%i.csv',folderName,macro_meso_iteration);
    csvwrite( outname,MacroExxColumn);
    
    % save the MacroEyyColumn
    %     outname = sprintf('./out%i/MacroEyyColumn%i.csv',folderNum,macro_meso_iteration);
    outname = sprintf('./%s/MacroEyyColumn%i.csv',folderName,macro_meso_iteration);
    csvwrite( outname,MacroEyyColumn);
    
    % save the MacroThetaColumn
    %     outname = sprintf('./out%i/MacroThetaColumn%i.csv',folderNum,macro_meso_iteration);
    outname = sprintf('./%s/MacroThetaColumn%i.csv',folderName,macro_meso_iteration);
    csvwrite( outname, MacroThetaColumn);
    
    % save the RhoColumn
    %     outname = sprintf('./out%i/RhoColumn%i.csv',folderNum,macro_meso_iteration);
    outname = sprintf('./%s/RhoColumn%i.csv',folderName,macro_meso_iteration);
    csvwrite( outname,  RhoColumn);
    %     end
    
    % subplot(xplots,yplots,c);c=c+1;
    % p.PlotArrayGeneric( ExxMacro, 'ExxMacro')
    %
    % subplot(xplots,yplots,c);c=c+1;
    % p.PlotArrayGeneric( EyyMacro, 'EyyMacro')
    %
    % subplot(xplots,yplots,c);c=c+1;
    % p.PlotArrayGeneric( ThetaMacro, 'ThetaMacro')
    %
    % x = [ExxMacro; EyyMacro; ThetaMacro]
    %
    % figure
    % subplot(2,2,1)
    % p.PlotArrayGeneric(100* relativeErrorExx, 'Percent Error Exx')
    % subplot(2,2,2)
    % p.PlotArrayGeneric( 100*relativeErrorEyy, 'Perecent Error Eyy')
    % subplot(2,2,3)
    % p.PlotArrayGeneric(100* relativeErrorTheta, 'Percent Error Theta')
    % subplot(2,2,4)
    % p.PlotArrayGeneric(diffTheta, 'Diff Theta')
  
    
    % --------------------------------------------------------
    %
    %    Compare Error in finding meso densities with their targets
    %    Meso Validation Case
    %
    % --------------------------------------------------------
  
    
    % --------------------------
    % Plot the raw data showing the density as circles
    % --------------------------
    RhoColor=RhoColumn; % color
    circleSize = ones(size(RhoColumn))*100; % circle size.
    scatter3(MacroExxColumn,MacroEyyColumn,MacroThetaColumn,circleSize,RhoColor,'filled','MarkerEdgeColor','k')
    title(sprintf('%s,Density Plot ',PlotType));
    colorbar
    xlabel('Exx');
    ylabel('Eyy');
    zlabel('Theta');
    %colormap('gray')
    % colormap(flipud(gray(256)));
    colormap('parula');
    
    nameGraph2 = sprintf('./%s_RawData%i.png',PlotType, config.macro_meso_iteration);
    print(nameGraph2,'-dpng');
    
    % ---------------------------
    % Calculate the error
    % ---------------------------
    errorType = 2; % relative
    if(errorType==1)
        % Exx Relative Error
        diffExx = MacroExxColumn-ExxArray;
        diffExx=abs(diffExx);
        diffExx=diffExx./MacroExxColumn; % Relative Error
        diffExx( MacroExxColumn==0)=0;
        
        % Eyy Relative Error
        diffEyy = MacroEyyColumn-EyyArray;
        diffEyy=abs(diffEyy);
        diffEyy=diffEyy./MacroEyyColumn; % Relative Error
        diffEyy( MacroEyyColumn==0)=0;
        
        % Theta Relative Error
        diffTheta = MacroThetaColumn-thetaArray;
        diffTheta=abs(diffTheta);
        diffTheta=diffTheta./MacroThetaColumn; % Relative Error
        diffTheta( MacroThetaColumn==0)=0;
        
        % calculate the total Error
        totalError = diffExx+(diffEyy)+(diffTheta);
        
        ErrorName = 'Relative';
        
    elseif(errorType==2)
        
        
        % Exx Relative DifferenceError
        diffExx = MacroExxColumn-ExxArray;
        diffExx=abs(diffExx);
        sumABSVAlues = abs(MacroExxColumn)+abs(ExxArray);
        diffExx=diffExx./(sumABSVAlues./2); %
        
%         size(MacroExxColumn)
%         size(diffExx)
        %         diffExx( MacroExxColumn==0)=0;
        
        % Eyy Relative Error
        diffEyy = MacroEyyColumn-EyyArray;
         diffEyy=abs(diffEyy);
        sumABSVAlues = abs(MacroEyyColumn)+abs(EyyArray);
        diffEyy=diffEyy./(sumABSVAlues./2); %
        %
        %         diffEyy( MacroEyyColumn==0)=0;
        
        % Theta Relative Error
        diffTheta = MacroThetaColumn-thetaArray;
         diffTheta=abs(diffTheta);
        sumABSVAlues = abs(MacroThetaColumn)+abs(thetaArray);
        diffTheta=diffTheta./(sumABSVAlues./2); %
        
        
        % calculate the total Error
        totalError = diffExx+(diffEyy)+(diffTheta);
        
        ErrorName = 'Relative Difference';
        
        
    end
    
    tttt=size(MacroExxColumn,1);
    if(tttt<1000)
        dotSize = 100;
    else
        dotSize = 20;
    end
    
    
    % --------------------------
    % Plot the Exx Error as circles
    % --------------------------
    figure
    
    ColorColumn=diffExx; % color
    circleSize = ones(size(ColorColumn))*dotSize; % circle size.
    scatter3(MacroExxColumn,MacroEyyColumn,MacroThetaColumn,circleSize,ColorColumn,'filled','MarkerEdgeColor','k')
    title(sprintf('Exx %s Error as circles (Target - Actual)/Target, %s',ErrorName,PlotType));
    colorbar
    xlabel('Exx');
    ylabel('Eyy');
    zlabel('Theta');
    caxis([0 1.5]);
    %colormap('gray')
    % colormap(flipud(gray(256)));
    colormap('parula');
    nameGraph2 = sprintf('./%s_ExxError%i.png', PlotType,config.macro_meso_iteration);
    print(nameGraph2,'-dpng');
    
    edges = [0:0.1:4];
    h = histogram(reshape(diffExx,1,[]),edges);
    title(sprintf('Exx %s Histogram %s',ErrorName,PlotType));
    nameGraph2 = sprintf('./%s_ExxError%iAsHistogram.png', PlotType,config.macro_meso_iteration);
    print(nameGraph2,'-dpng');
    
    ExxSummedError = sum(sum(diffExx));
    ExxAvgError = mean(reshape(diffExx,1,[]));
    ExxSTDError = std(reshape(diffExx,1,[]));
    ExxMedian = median(reshape(diffExx,1,[]));
    ExxMode = mode(reshape(diffExx,1,[]));
    fprintf('\n--------\nSummed Exx Error %f Average Error %f and STD %f Median  %f, and Mode %f\n----\n',ExxSummedError,ExxAvgError,ExxSTDError,ExxMedian,ExxMode);
    
    % --------------------------
    % Plot the Eyy Error as circles
    % --------------------------
    figure
    
    ColorColumn=diffEyy; % color
    circleSize = ones(size(ColorColumn))*dotSize; % circle size.
    scatter3(MacroExxColumn,MacroEyyColumn,MacroThetaColumn,circleSize,ColorColumn,'filled','MarkerEdgeColor','k')
    title(sprintf('Eyy %s Error as circles (Target - Actual)/Target,%s',ErrorName,PlotType));
    colorbar
    xlabel('Exx');
    ylabel('Eyy');
    zlabel('Theta');
    %colormap('gray')
    % colormap(flipud(gray(256)));
    colormap('parula');
    caxis([0 1.5]);
    nameGraph2 = sprintf('./%s_EyyError%i.png', PlotType,config.macro_meso_iteration);
    print(nameGraph2,'-dpng');
    EyySummedError = sum(sum(diffEyy));
    EyyAvgError = mean(reshape(diffEyy,1,[]));
    EyySTDError = std(reshape(diffEyy,1,[]));
    EyyMedian = median(reshape(diffEyy,1,[]));
    EyyMode = mode(reshape(diffEyy,1,[]));
    fprintf('\n--------\nSummed Eyy Error %f Average Error %f and STD %f Median %f, and Mode %f\n----\n',EyySummedError,EyyAvgError,EyySTDError,EyyMedian,EyyMode);
    
    
    h = histogram(reshape(diffEyy,1,[]),edges);
    title(sprintf('Eyy %s Error Histogram,%s',ErrorName,PlotType));
    nameGraph2 = sprintf('./%s_EyyError%iAsHistogram.png',PlotType, config.macro_meso_iteration);
    print(nameGraph2,'-dpng');
    
    
    % --------------------------
    % Plot the Theta Error as circles
    % --------------------------
    figure
    
    ColorColumn=diffTheta; % color
    circleSize = ones(size(ColorColumn))*dotSize; % circle size.
    scatter3(MacroExxColumn,MacroEyyColumn,MacroThetaColumn,circleSize,ColorColumn,'filled','MarkerEdgeColor','k')
    title(sprintf('Theta %s Error as circles (Target - Actual)/Target, %s',ErrorName,PlotType));
    colorbar
    xlabel('Exx');
    ylabel('Eyy');
    zlabel('Theta');
    caxis([0 1.5]);
    %colormap('gray')
    % colormap(flipud(gray(256)));
    colormap('parula');
    nameGraph2 = sprintf('./%s_ThetaError%i.png', PlotType,config.macro_meso_iteration);
    print(nameGraph2,'-dpng');
    ThetaSummedError = sum(sum(diffTheta));
    ThetaAvgError = mean(reshape(diffTheta,1,[]));
    ThetaSTDError = std(reshape(diffTheta,1,[]));
    ThetaMedian = median(reshape(diffTheta,1,[]));
    ThetaMode = mode(reshape(diffTheta,1,[]));
    fprintf('\n--------\nSummed Theta Error %f Average Error %f and STD %f Median %f and Mode %f \n----\n',ThetaSummedError,ThetaAvgError,ThetaSTDError,ThetaMedian,ThetaMode);
    
    
    
    h = histogram(reshape(diffTheta,1,[]),edges);
    title(sprintf('Theta %s Error Histogram, %s',ErrorName,PlotType));
    nameGraph2 = sprintf('./%s_ThetaError%iAsHistogram.png',PlotType, config.macro_meso_iteration);
    print(nameGraph2,'-dpng');
    
    
    % --------------------------
    % Plot the Combined Normalized Error
    % --------------------------
    figure
    % take ABS, and normalize
    %         totalError = abs(diffExx)/matProp.E_material1+abs(diffEyy)/matProp.E_material1+abs(diffTheta)/(pi/2);
    
    ColorColumn=totalError; % color
    circleSize = ones(size(ColorColumn))*dotSize; % circle size.
    scatter3(MacroExxColumn,MacroEyyColumn,MacroThetaColumn,circleSize,ColorColumn,'filled','MarkerEdgeColor','k')
    title(sprintf('Summed Error as circles for Exx, Eyy, Theta, %s',PlotType));
    colorbar
    xlabel('Exx');
    ylabel('Eyy');
    zlabel('Theta');
    caxis([0 3]);
    %colormap('gray')
    % colormap(flipud(gray(256)));
    colormap('parula');
    nameGraph2 = sprintf('./%s_combinedError%i.png', PlotType,config.macro_meso_iteration);
    print(nameGraph2,'-dpng');
    
    
    
    TotalSummedError = sum(sum(totalError));
    TotalAvgError = mean(reshape(totalError,1,[]));
    TotalSTDError = std(reshape(totalError,1,[]));
    TotalErrorMedian = median(reshape(totalError,1,[]));
    TotalErrorMode = mode(reshape(totalError,1,[]));
    fprintf('\n--------\nSummed totalError Error %f Average Error %f and STD %f, Median %f, and Mode %f\n----\n',TotalSummedError,TotalAvgError,TotalSTDError,TotalErrorMedian,TotalErrorMode);
    
    h = histogram(reshape(totalError,1,[]),edges);
    title(sprintf(' Summed Error Histogram, %s',PlotType));
    nameGraph2 = sprintf('./%son_combinedHistorgramError%iAsHistogram.png',PlotType, config.macro_meso_iteration);
    print(nameGraph2,'-dpng');
    
    
    statData = [ExxSummedError;ExxAvgError;ExxSTDError;ExxMedian;ExxMode; ...
        EyySummedError;EyyAvgError;EyySTDError;EyyMedian;EyyMode; ...
        ThetaSummedError;ThetaAvgError;ThetaSTDError;ThetaMedian;ThetaMode; ...
        TotalSummedError;TotalAvgError;TotalSTDError;TotalErrorMedian;TotalErrorMode]
    csvwrite(sprintf('./%s_Metrics.csv',PlotType),statData);
    
    fprintf('MesoMethod %i (0=feedback, 1 lookup) and Material Update Scheme %i\n',config.UseLookUpTableForPsuedoStrain,config.mesoVolumeUpdateMethod);
    fprintf('Lookup Search Method %i (2 = search table, scale eta, 4 = particle swarmn\n',config.lookupSearchScheme);
    
    % Change locations where we the target was a zero, to be an
    % average error value. The error is not acutally zero, so we
    % dont' want to imply this.
    
    diffExx( MacroExxColumn==0)=ExxAvgError;
    diffEyy( MacroEyyColumn==0)=EyyAvgError;
    diffTheta( MacroThetaColumn==0)=ThetaAvgError;
    totalErrorV2=diffExx+diffEyy+diffTheta;
    %     totalErrorV2=reshape(totalErrorV2,[],1);
    %     totalErrorV2=totalErrorV2';
    %     size(MacroExxColumn)
    %     size(totalErrorV2)
    
    
    
    if(config.UseLookUpTableForPsuedoStrain==1)
        m = config.mesoVolumeUpdateMethod;
    else
        m=config.lookupSearchScheme;
    end
    outname = sprintf('./%s/%s_Valid_TotalError_Method%i_Config%i.csv',folderName,PlotType,config.UseLookUpTableForPsuedoStrain,m);
    csvwrite(outname,totalErrorV2);
    
    %     end
    
    % ----------------------
    % ---------------------
    
    % Plot the metrics over several iterations
    if(config.macro_meso_iteration>1 || config.validationModeOn==1)
        AllData=[];
        MacroExxColumnTotal=[];
        MacroEyyColumnTotal=[];
        MacroThetaColumnTotal=[];
        RhoColumnTotal=[];
        ErrorTotal=[];
        for jj = 1:config.macro_meso_iteration
            if(config.validationModeOn==1)
                PlotType = 'MesoValidation' ;
            else
                PlotType=sprintf('Results From Iter %i',jj);
            end
            %              outname = sprintf('./%s/%s_Valid_TotalError_Method%i_Config%i.csv',folderName,PlotType,config.UseLookUpTableForPsuedoStrain,m);
            
            statData=  csvread(sprintf('./%s_Metrics.csv',PlotType));
            
            
            AllData=[AllData statData];
            
            
            
            % ----------
            % ----------
            macro_meso_iteration=jj;
            outnameExx = sprintf('./%s/MacroExxColumn%i.csv',folderName,macro_meso_iteration);
            %             outname = sprintf('./%s/MacroExxColumn%i.csv',folderName,macro_meso_iteration);
            outnameEyy = sprintf('./%s/MacroEyyColumn%i.csv',folderName,macro_meso_iteration);
            outnametheta = sprintf('./%s/MacroThetaColumn%i.csv',folderName,macro_meso_iteration);
            outnamerho = sprintf('./%s/RhoColumn%i.csv',folderName,macro_meso_iteration);
            outnameTotalError = sprintf('./%s/%s_Valid_TotalError_Method%i_Config%i.csv',folderName,PlotType,config.UseLookUpTableForPsuedoStrain,m);
            
            MacroExxColumn2=csvread(outnameExx);
            MacroExxColumnTotal=[MacroExxColumnTotal; MacroExxColumn2];
            %             fprintf('size of MacroExxColumn data %i\n',jj);
            %             size(MacroExxColumn2)
            
            % save the MacroEyyColumn
            temp=csvread(outnameEyy);
            MacroEyyColumnTotal=[MacroEyyColumnTotal; temp];
            
            % save the MacroThetaColumn%
            temp=csvread(outnametheta);
            MacroThetaColumnTotal=[MacroThetaColumnTotal; temp];
            
            % save the RhoColumn
            temp=csvread(outnamerho);
            RhoColumnTotal=[RhoColumnTotal; temp];
            
            temp=csvread(outnameTotalError);
            %             fprintf('size of error data %i\n',jj);
            %             size(temp)
            ErrorTotal=[ErrorTotal; temp];
            
            
        end
        
        x=1:size(AllData,2);
        
        for kk = 1:20
            plot(x,AllData(kk,:))
            hold on
        end
        %         AllData
        csvwrite(sprintf('%s_allData.csv',PlotType),AllData);
        
        
        legend('ExxSummedError','ExxAvgError','ExxSTDError','ExxMedian','ExxMode', ...
            ' EyySummedError','EyyAvgError','EyySTDError','EyyMedian','EyyMode', ...
            ' ThetaSummedError','ThetaAvgError','ThetaSTDError','ThetaMedian','ThetaMode', ...
            'TotalSummedError','TotalAvgError','TotalSTDError','TotalErrorMedian','TotalErrorMode');
        nameGraph2 = sprintf('./%s_MetricsPlotOverSevearlIterations.png',PlotType);
        print(nameGraph2,'-dpng');
        hold off
        
        figure(2)
        
        
        
        %         size(MacroExxColumnTotal)
        %         size(MacroEyyColumnTotal)
        %         size(MacroThetaColumnTotal)
        %         size(ErrorTotal)
    
        
        
        % Save the combined Erorr Values
        folderName='ErrorData';
        outname = sprintf('./%s/MacroExxColumnTotal%i.csv',folderName,macro_meso_iteration);
        csvwrite( outname,MacroExxColumnTotal);
        
        % save the MacroEyyColumn
        %     outname = sprintf('./out%i/MacroEyyColumn%i.csv',folderNum,macro_meso_iteration);
        outname = sprintf('./%s/MacroEyyColumnTotal%i.csv',folderName,macro_meso_iteration);
        csvwrite( outname,MacroEyyColumnTotal);
        
        % save the MacroThetaColumn
        %     outname = sprintf('./out%i/MacroThetaColumn%i.csv',folderNum,macro_meso_iteration);
        outname = sprintf('./%s/MacroThetaColumnTotal%i.csv',folderName,macro_meso_iteration);
        csvwrite( outname, MacroThetaColumnTotal);
        
        % save the RhoColumn
        %     outname = sprintf('./out%i/RhoColumn%i.csv',folderNum,macro_meso_iteration);
        outname = sprintf('./%s/ErrorTotal%i.csv',folderName,macro_meso_iteration);
        csvwrite( outname,  ErrorTotal);
        
            % -----------------------------
        % -----------------------------
        % -----------------------------
        % Plot the Total error
        % -----------------------------
        % -----------------------------
        % -----------------------------
        
        errorIncrements= [2];
        
        azArray=[  0  90  0 ];
        elArray=[  90 0   180 ];
        for eIncrement = errorIncrements
            fprintf('plotting values from all iterations with error more than %i\n',eIncrement*100);
            logic = ErrorTotal>eIncrement;
            MacroExxColumnTotal=MacroExxColumnTotal(logic);
            MacroEyyColumnTotal=MacroEyyColumnTotal(logic);
            MacroThetaColumnTotal=MacroThetaColumnTotal(logic);
            ErrorTotal=ErrorTotal(logic);
            
            for kkkkk = 1:3
                az=azArray(kkkkk);
                el = elArray(kkkkk);
                subplot(2,2,kkkkk)
                ColorColumn=ErrorTotal; % color
                circleSize = ones(size(ColorColumn))*10; % circle size.
                scatter3(MacroExxColumnTotal,MacroEyyColumnTotal,MacroThetaColumnTotal,circleSize,ColorColumn,'filled')
                title(sprintf(' Summed Error %s greaterThan %i',PlotType,eIncrement*100));
                colorbar
                xlabel('Exx');
                ylabel('Eyy');
                zlabel('Theta');
                caxis([0 3]);
                %colormap('gray')
                % colormap(flipud(gray(256)));
                colormap('parula');
                view(az, el);
                
            end
            
            nameGraph2 = sprintf('./%s_combinedError%i_AllData4by4_largethan%i.png', PlotType,config.macro_meso_iteration,eIncrement*100);
            print(nameGraph2,'-dpng', '-r600')
            
            subplot(1,1,1)
            ColorColumn=ErrorTotal; % color
            circleSize = ones(size(ColorColumn))*10; % circle size.
            scatter3(MacroExxColumnTotal,MacroEyyColumnTotal,MacroThetaColumnTotal,circleSize,ColorColumn,'filled')
            title(sprintf(' Summed Error %s greater than %i',PlotType,eIncrement*100));
            colorbar
            xlabel('Exx');
            ylabel('Eyy');
            zlabel('Theta');
            caxis([0 3]);
            %colormap('gray')
            % colormap(flipud(gray(256)));
            colormap('parula');
            nameGraph2 = sprintf('./%s_combinedError%i_AllDataJustOne_largerthan_%i.png', PlotType,config.macro_meso_iteration,eIncrement*100);
            print(nameGraph2,'-dpng', '-r600')
        end
        
        
    end
end




% -----------------------------------
%
%    Generate surface fit.
%   must be commented out for matlab to compile on cluster.
%
% % -----------------------------------
% if(1==1)
%     MacroExxColumnTotal=[];
%     MacroEyyColumnTotal=[];
%     MacroThetaColumnTotal=[];
%     RhoColumnTotal=[];
%     folderName='data';
%
%     useSubSysValues = 1;
%
%     if(useSubSysValues==0)
%         outnameExx = sprintf('./%s/MacroExxColumn%i.csv',folderName,macro_meso_iteration);
%         outnameEyy = sprintf('./%s/MacroEyyColumn%i.csv',folderName,macro_meso_iteration);
%         outnametheta = sprintf('./%s/MacroThetaColumn%i.csv',folderName,macro_meso_iteration);
%         outnamerho = sprintf('./%s/RhoColumn%i.csv',folderName,macro_meso_iteration);
%     else
%         outnameExx = sprintf('./%s/ExxSubSysValues%i.csv',folderName,macro_meso_iteration);
%         outnameEyy = sprintf('./%s/EyySubSysValues%i.csv',folderName,macro_meso_iteration);
%         outnametheta = sprintf('./%s/ThetaSubSysValues%i.csv',folderName,macro_meso_iteration);
% %         outnamerho = sprintf('./%s/densityUsedSubSysValues%i.csv',folderName,macro_meso_iteration);
%          outnamerho = sprintf('./%s/RhoColumn%i.csv',folderName,macro_meso_iteration);
%     end
%
%     for i = 1:macro_meso_iteration
%         %------------------------
%         % read the macro columns as well. This will help with future analysis
%         % ----------------------------
%         % save the MacroExxColumn
%         MacroExxColumn=csvread(outnameExx);
%         MacroExxColumnTotal=[MacroExxColumnTotal; MacroExxColumn];
%
%         % save the MacroEyyColumn
%         temp=csvread(outnameEyy);
%         MacroEyyColumnTotal=[MacroEyyColumnTotal; temp];
%
%         % save the MacroThetaColumn%
%         temp=csvread(outnametheta);
%         MacroThetaColumnTotal=[MacroThetaColumnTotal; temp];
%
%         % save the RhoColumn
%         temp=csvread(outnamerho);
%         RhoColumnTotal=[RhoColumnTotal; temp];
%     end
%
% %     for jjj=1:5
% %     % Add full dense case
% %       MacroExxColumnTotal=[MacroExxColumnTotal; max(MacroExxColumnTotal)];
% %       MacroEyyColumnTotal=[MacroEyyColumnTotal; max(MacroEyyColumnTotal)];
% %       MacroThetaColumnTotal=[MacroThetaColumnTotal; 0];
% %       RhoColumnTotal=[RhoColumnTotal;1];
% %     end
%
%
%     x=MacroExxColumnTotal/matProp.E_material1;
%     y = MacroEyyColumnTotal/matProp.E_material1;
%
%      if(useSubSysValues==1)
%         x=x';
%         y = y';
%         MacroThetaColumnTotal=MacroThetaColumnTotal';
%      end
%
%     z=RhoColumnTotal;
%
%    options= fitoptions;
% %    options.Normalize ='on';
% %    options.fittype='poly22';
%      f1 = fit([x y],z,'poly33',options)
%      coeffvalues(f1)
% %      f2 = fit([x y],z,'poly23', 'Exclude', z > 1);
% %     o=Optimizer;
% %     [~, ~,annZ] = o.CalculateDensitySensitivityandRho(x,y,MacroThetaColumnTotal,ones(size(MacroEyyColumnTotal)),DV.ResponseSurfaceCoefficents,config,matProp,0);
%
%     figure
%      plot(f1, [x y], z);
% %     hold on
% %     scatter3(x,y,z,'b')
% %     hold on
% %     scatter3(x,y,annZ,'r');
%
%     title('Fit with data points. Red=Ann, Blue=Actual ')
%     xlabel('Exx');
%     ylabel('Eyy');
%     zlabel('rho');
%     zlim([0 1])
%     size(RhoColumnTotal)
% end


%annTest(macro_meso_iteration);

% if 1==0
%     Exx = matProp.E_material1;
%     Eyy = 0;
%     v= 1;
%     rhoArray = [rhoArray;v];
%     ExxArray=[ExxArray;Exx];
%     EyyArray=[EyyArray;Eyy];
%
%      Exx = matProp.E_material1;dos
%     Eyy = matProp.E_material1/2;
%     v= 1;
%     rhoArray = [rhoArray;v];
%     ExxArray=[ExxArray;Exx];
%     EyyArray=[EyyArray;Eyy];
%
%     % Both extremes
%        Exx = matProp.E_material1;
%     Eyy = matProp.E_material1;
%     v= 1;
%     rhoArray = [rhoArray;v];
%     ExxArray=[ExxArray;Exx];
%     EyyArray=[EyyArray;Eyy];
%
%       Exx = matProp.E_material1/2;
%     Eyy =matProp.E_material1 ;
%     v= 1;
%     rhoArray = [rhoArray;v];
%     ExxArray=[ExxArray;Exx];
%     EyyArray=[EyyArray;Eyy];
%
%     %Eyy extreme
%     Exx =0;
%     Eyy =  matProp.E_material1;
%     v= 1;
%     rhoArray = [rhoArray;v];
%     ExxArray=[ExxArray;Exx];
%     EyyArray=[EyyArray;Eyy];
%
%
% end


% if(1==1)
%
%     nameArray = sprintf('./out%i/ExxArrayForFitting%i.csv',folderNum, config.macro_meso_iteration);
%     csvwrite(nameArray,ExxArray);
%
%     nameArray = sprintf('./out%i/EyyArrayForFitting%i.csv',folderNum, config.macro_meso_iteration);
%     csvwrite(nameArray,EyyArray);
%
%     nameArray = sprintf('./out%i/ThetaArrayForFitting%i.csv',folderNum, config.macro_meso_iteration);
%     csvwrite(nameArray,thetaArray);
%
%
%
%     nameArray = sprintf('./out%i/RhoArrayForFitting%i.csv',folderNum, config.macro_meso_iteration);
%     csvwrite(nameArray,rhoArray);
%
%     % REad the old arrays as well.
%     if(config.macro_meso_iteration>1)
%         for jjj= 1:config.macro_meso_iteration-1
%             nameArray = sprintf('./out%i/ExxArrayForFitting%i.csv',folderNum, jjj);
%             MacroExxColumnTemp =  csvread(nameArray);
%             ExxArray=[ExxArray ;MacroExxColumnTemp];
%
%             nameArray = sprintf('./out%i/EyyArrayForFitting%i.csv',folderNum, jjj);
%             MacroEyyColumnTemp =  csvread(nameArray);
%             EyyArray=[EyyArray; MacroEyyColumnTemp];
%
%
%             nameArray = sprintf('./out%i/ThetaArrayForFitting%i.csv',folderNum, jjj);
%             thetaArrayTemp =  csvread(nameArray);
%             thetaArray=[thetaArray; thetaArrayTemp];
%
%
%
%             nameArray = sprintf('./out%i/RhoArrayForFitting%i.csv',folderNum, jjj);
%             rhoArrayTemp =  csvread(nameArray);
%             rhoArray=[rhoArray; rhoArrayTemp];
%         end
%
%
%     end
%
%
%
%
%     figure(1)
%     scaleUp = matProp.E_material1;
%     config.useThetaInSurfaceFit=1;
%     if(config.useThetaInSurfaceFit==1)
%
%         % -----------------------
%         % PLot the raw data. Scale up the rho, to make the color range
%         % larger.
%         % -----------------------
%         RhoColor=rhoArray; % color
%         circleSize = ones(size(ExxArray))*100; % circle size.
%         scatter3(ExxArray,EyyArray,thetaArray,circleSize,RhoColor);
%         xlabel('Exx');
%         ylabel('Eyy');
%         zlabel('Theta');
%         title(sprintf('Rho (the color) as a function of Exx, Eyy, theta: iter %i',config.macro_meso_iteration));
%         colorbar
%
%         nameGraph2 = sprintf('./RhoDensityOfExxEyyPlot%i.png', config.macro_meso_iteration);
%         print(nameGraph2,'-dpng');
%
%         % ----------------------------
%         % Plot 2, with symmmetry taken into account.
%         % Make Exx always larger
%         % theta between 0 and pi/4
%         % ----------------------------
%
%       % Make the inputs be so taht Exx > Eyy
% % Rather than a strict theta, use the distance from pi/4, since the problem
% % is symmetric arround pi/4
% temp = ExxArray;
% logic = EyyArray>ExxArray;
% ExxArray(logic)=EyyArray(logic);
% EyyArray(logic) =temp(logic);
% %
% % min(thetaArray)
% % max(thetaArray)
% % thetaArray=((pi/4)^2+thetaArray.^2).^(1/2);
% temp2 = thetaArray;
% logic2 = thetaArray>pi/4;
% logic3 = thetaArray<pi/4;
% thetaArray(logic2)=thetaArray(logic2)-pi/4;
% thetaArray(logic3)=pi/4-thetaArray(logic3);
%
%          figure
%           scatter3(ExxArray,EyyArray,thetaArray,circleSize,RhoColor);
%         xlabel('Exx');
%         ylabel('Eyy');
%         zlabel('Theta');
%         title(sprintf('Plot 2Rho (the color) as a function of Exx, Eyy, theta: iter %i',config.macro_meso_iteration));
%         colorbar
%
%           nameGraph2 = sprintf('./Plot2RhoDensityOfExxEyyPlot%i.png', config.macro_meso_iteration);
%         print(nameGraph2,'-dpng');
%
%
%         % -------------------------------
%         % Least squares fit
%         % ------------------------------
% %
% x0=ones(1,10);
% x0 = randi([-5,5],1,10);
% A = [];
% b = [];
%
%
% % theta the same
% % rho the same.
% % scale down Exx, Eyy
%
% X = ExxArray/scaleUp;
% Y = EyyArray/scaleUp;
%
% %         Z = thetaArray/(pi/4);
% Z = thetaArray;
% R = rhoArray;
%
%
% ub = ones(6,1)*10000;
% lb = -ub;
% o=Optimizer;
% [coefficients finalObjective]= fmincon(@(x) fitObjectiveV2(x,X,Y,Z,R,o,config,matProp),x0,A,b,[],[],lb,ub);
%
% finalObjective
%
% % use the scaled data
% %         numPointsXandY = 20;
% %tt  =1/numPointsXandY;
% tt  =0.05;
%
% [Xgrid, Ygrid, Zgrid]=meshgrid(0:tt:max(X),0:tt:max(Y),0:0.2:max(Z));
%
% % Reshape into columns
% E_xx=reshape(Xgrid,[],1);
% E_yy=reshape(Ygrid,[],1);
% theta=reshape(Zgrid,[],1);
%
% x=coefficients;
%
% %------------
% % Calcualte the rho values using the fitting polynomial
% %------------
% [~, ~,rhoExperimental] = o.CalculateDensitySensitivityandRho(E_xx,E_yy,theta,coefficients,config,matProp);

% -----------------------
% Plot
% - rescale the data to the correct form
% - plot scatter3
%-------------------------
%         rhoExperimental=rhoExperimental/scaleUp;
%         E_xx=E_xx*scaleUp;
%         E_yy=E_yy*scaleUp;
%    rhoExperimental=rhoExperimental;
%    theta=theta;


% %         E_xx(E_yy>E_xx)=0;
% %         E_yy(E_yy>E_xx)=0;
%
%
%         figure
%         circleSize = ones(size(E_xx))*100; % circle size.
%         scatter3(E_xx,E_yy,theta,circleSize,rhoExperimental,'filled');
%         title(sprintf('Response surface,Rho (the color) as a function of Exx, Eyy, theta: iter %i',config.macro_meso_iteration));
%         colorbar
%          xlabel('Exx');
%         ylabel('Eyy');
%         zlabel('Theta');
%         %                 hold off
%
%         nameGraph2 = sprintf('./RhoDensityOfExxEyyThetaResponseSurfacePlot%i.png', config.macro_meso_iteration);
%         print(nameGraph2,'-dpng');
%
%           nameArray = sprintf('./out%i/ExxEyyRhoFitCoefficients%i.csv',folderNum, config.macro_meso_iteration);
%         %      csvwrite(nameArray,sfArray);
%
%         dlmwrite(nameArray, x, 'delimiter', ',', 'precision', 15);
%
%
%     else
%         % ---------------------------------------------------
%         % Plot and surface fit where
%         %
%         % rho is a function of Exx and Eyy
%         % ---------------------------------------------------
%         scatter3(ExxArray,EyyArray,rhoArray);
%         xlabel('Exx');
%         ylabel('Eyy');
%         zlabel('Rho,Density');
%
%         x0=[1 1 1 1 1 1];
%         A = [];
%         b = [];
%
%         scaleUp = matProp.E_material1;
%
%         %      X = MacroExxColumn/matProp.E_material1;
%         %      Y = MacroEyyColumn/matProp.E_material1;
%         X = ExxArray;
%         Y = EyyArray;
%
%         Z = rhoArray*scaleUp;
%
%         ub = ones(6,1)*100000;
%         lb = -ub;
%         coefficients= fmincon(@(x) fitObjective(x,X,Y,Z),x0,A,b,[],[],lb,ub);
%         %       [x,fval,exitflag] = ga(@(x) fitObjective(x,X,Y,Z),  6,A,b,[],[],lb,ub);
%         %        sfArray=x;
%         coefficients=coefficients/scaleUp
%         %
%         [Xgrid, Ygrid]=meshgrid(0:1000:max(X),0:1000:max(Y));
%
%         %      Xgrid=Xgrid*matProp.E_material1;
%         %      Ygrid=Ygrid*matProp.E_material1;
%
%         p00 =coefficients(1);
%         p10=coefficients(2);
%         p01=coefficients(3) ;
%         p20=coefficients(4);
%         p11=coefficients(5);
%         p02=coefficients(6);
%
%         Zexperimental=   p00 + p10*Xgrid + p01*Ygrid + p20*Xgrid.^2 + p11*Xgrid.*Ygrid + p02*Ygrid.^2;
%
%         hold on
%         surf(Xgrid,Ygrid,Zexperimental)
%         hold off
%
%         nameArray = sprintf('./out%i/ExxEyyRhoFitCoefficients%i.csv',folderNum, config.macro_meso_iteration);
%         %      csvwrite(nameArray,sfArray);
%
%         dlmwrite(nameArray, coefficients, 'delimiter', ',', 'precision', 15);
%
%
%         nameGraph2 = sprintf('./RhoDensityOfExxEyyPlot%i.png', config.macro_meso_iteration);
%         print(nameGraph2,'-dpng');
%     end
%
%
%
%
%
%
%
%
%
%     %           sf = fit([MacroExxColumn, MacroEyyColumn],rhoArray,'poly22')
%     %           plot(sf,[MacroExxColumn,MacroEyyColumn],rhoArray)
%     %          xlabel('MacroExxColumn');
%     %     ylabel('MacroEyyColumn');
%     %     zlabel('Rho,Density');
%     %      sfArray = [sf.p00 sf.p10 sf.p01 sf.p20 sf.p02 sf.p11 ];
%
%
%
%
%
%     %       figure(3)
%     %      sf = fit([ExxArray, EyyArray],rhoArray,'poly23')
%     %      plot(sf,[ExxArray,EyyArray],rhoArray)
%     %         xlabel('ExxArray');
%     %     ylabel('EyyArray');
%     %     zlabel('Rho,Density');
%     % %
%     %
%     %     xlabel('Exx');
%     %     ylabel('Eyy');
%     %     zlabel('Rho,Density');
%     %
%     %     figure(3)
%     %     scatter3(ExxArray,EyyArray,thetaArray);
%     %     xlabel('Exx');
%     %     ylabel('Eyy');
%     %     zlabel('theta');
%
%     %     figure(4)
%     %        scatter3(MacroExxColumn,MacroEyyColumn,rhoArray);
%     %     histogram(thetaArray)
%
%     %     figure(5)
%     %     DV = DV.CalculateVolumeFractions( config,matProp);
% else
%     DV.CalculateVolumeFractions( config,matProp)
%     p = plotResults;
%     FEACalls=1;
%     p.plotTopAndFraction(DV,  config, matProp, FEACalls); % plot the results.
%
%
%     nameGraph = sprintf('./MesoDesignExxEyyThetaVars%i.png', config.macro_meso_iteration);
%     print(nameGraph,'-dpng');
% end
% %

% validationMeso =1;
% if(validationMeso ==1)
%
%     totalValidationProblems=config.nely*config.nelx;
%     numSegments = floor(totalValidationProblems^(1/3));
%     numSegmentsExx = numSegments;
%     numSegmentsTheta = floor(totalValidationProblems/(numSegmentsExx^2));
%
%     numSegmentsExx=numSegmentsExx-1;
%     numSegmentsTheta=numSegmentsTheta-1;
%
%     ExxVector =0:matProp.E_material1/numSegmentsExx:matProp.E_material1;
%     EyyVector =0:matProp.E_material1/numSegmentsExx:matProp.E_material1;
%     thetaVector = 0:(pi/2)/numSegmentsTheta:pi/2;
%     [ExxValues EyyValues ThetaValues] = meshgrid(ExxVector,ExxVector,thetaVector);
%
%     %ExxValues=padarray(ExxValues,
%     ExxValues=reshape(ExxValues,1,[]);
%     [t1 t2]=size(ExxValues);
%     ExxValues=padarray(ExxValues,[0 totalValidationProblems-t2],'post');
%     ExxValues=reshape(ExxValues,config.nely,config.nelx);
%
%     EyyValues=reshape(EyyValues,1,[]);
%     EyyValues=padarray(EyyValues,[0 totalValidationProblems-t2],'post');
%     EyyValues=reshape(EyyValues,config.nely,config.nelx);
%
%     ThetaValues=reshape(ThetaValues,1,[]);
%     ThetaValues=padarray(ThetaValues,[0 totalValidationProblems-t2],'post');
%     ThetaValues=reshape(ThetaValues,config.nely,config.nelx);
%
%     xSimp = ones(t1, t2);
%     xSimp=reshape(xSimp,1,[]);
%     xSimp=padarray(xSimp,[0 totalValidationProblems-t2],'post');
%     xSimp=reshape(xSimp,config.nely,config.nelx);
%
%
%     ExxNewArray = [];
%      EyyNewArray = [];
%       EyyNewArray = [];
%        EzzNewArray = [];
%     for e = 1:t1*t2
%         Xcondition = ExxVector(1)==DV.Exx(e);
%         Ycondition = 1;%ExxVector(1)==DV.Exx(e);
%         thetaCondtion =1;% ThetaValues(1)==DV.t(e);
%         %rhoCondtion = ExxVector(1)==DV.Exx(e);
%         if(Xcondition==1 && Ycondition==1 && thetaCondtion==1)
%
%         end
%
%     end
%
%
% end



