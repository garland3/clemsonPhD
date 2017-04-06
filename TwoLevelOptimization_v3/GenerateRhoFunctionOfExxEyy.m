function []=GenerateRhoFunctionOfExxEyy(config)
postProcess = 1;
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



% save the Exx field
outname = sprintf('./out%i/ExxValues%i.csv',folderNum,macro_meso_iteration);
ExxMacro = csvread(outname);

% save the Eyy field
outname = sprintf('./out%i/EyyValues%i.csv',folderNum,macro_meso_iteration);
EyyMacro =csvread(outname);

% save the Theta field
outname = sprintf('./out%i/ThetaValues%i.csv',folderNum,macro_meso_iteration);
ThetaMacro = csvread(outname);



matProp=MaterialProperties;
DV.x = xxx;
rhoArray = [];
ExxArray = [];
EyyArray = [];
thetaArray=[];

for e = 1:ne %ne:-1:1
    fprintf('element %i of %i\n',e,ne);
    macroElementProps.elementNumber=e;
    results = elementXYposition(macroElementProps.elementNumber,:);
    macroElementProps.yPos = results(1);
    macroElementProps.xPos = results(2);
    elementNumber=e;
    macroElementProps.densitySIMP = xxx(macroElementProps.yPos,macroElementProps.xPos );
    ActualThetaValue = ThetaMacro(macroElementProps.yPos,macroElementProps.xPos );
    ActualExx = ExxMacro(macroElementProps.yPos,macroElementProps.xPos );
    ActualEyy = EyyMacro(macroElementProps.yPos,macroElementProps.xPos );
    
  
    
    
    if(macroElementProps.densitySIMP>config.voidMaterialDensityCutOff)
        % save the psuedo strain values
        %         outname = sprintf('./out%i/psuedostrain_Ite%i_forElement%i.csv',folderNum,macro_meso_iteration,elementNumber);
        %     p =  macroElementProperties.psuedoStrain;
        %               p =  csvread(outname);
        
        % save the final volume
        outname = sprintf('./out%i/volumeUsed_Ite%i_forElement%i.csv',folderNum,macro_meso_iteration,elementNumber);
        %              v =    configMeso.totalVolume;
        %     csvwrite(outname,v);
        if exist(outname, 'file') ~= 2
            continue;
        end
        v =  csvread(outname);
        
        
        outname = sprintf('./out%i/Dmatrix_%i_forElement_%i.csv',folderNum,macro_meso_iteration,elementNumber);
%          outname = sprintf('./out%i/DsystemIter%i_Element_%i.csv',folderNum,macro_meso_iteration,elementNumber);
        Din = csvread(outname);
        
            Dcalculated= matProp.getDmatMatrixTopExxYyyRotVars(config, macroElementProps.densitySIMP,ActualExx, ActualEyy,ActualThetaValue,1);
       
    
%         temp33=pi/2-ActualThetaValue;
%         incr = pi/360;
%         thetaValuesToTest1=temp33-4*pi/8+incr:incr:temp33+4*pi/8-incr;
          
     
        
   
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
             
       
        
        % in the case  where, Exx = Eyy, then the material is basically
        % isotropic and rotating to find the orthotropic orientaiton will
        % not work. In this case, set the theta to the actualTheta        
        % the criteria is that Exx and Eyy must be within 0.2% of each
        % other's value.
        if(abs(100*(Exx-Eyy)/Exx)<0.2)
            Theta=ActualThetaValue;
            str = sprintf('Exx = Eyy, setting theta to sys theta')
        end
        
        % Also, there is some difficulty when actualTheta is 0 or pi/2
        diffTheta = abs( ActualThetaValue-Theta);
        
        if(diffTheta>(pi/2-epsilon))
            if(Theta>pi/4)
                Theta=pi/2-Theta;
            else
                Theta=pi/2-Theta;
            end
            
              str = sprintf('Theta on wrong boundary. Switching values. ')
                diffTheta = abs( ActualThetaValue-Theta); % The new Diff Theta
        end
   
        
        
         diffX = ActualExx-Exx;
        diffY = ActualEyy-Eyy;
        
        
        relativeErrorDiffExx = diffX/ActualExx;
         relativeErrorDiffyy = diffY/ActualEyy;
         relativeErrorDiffTheta = diffTheta/ActualThetaValue;
        if(relativeErrorDiffExx>0.4 || relativeErrorDiffyy>0.4  || relativeErrorDiffTheta>0.4)
            sprintf('%i Large Error', e)
        end
    
        
        
        
        
        %DV.x already saved.
        DV.Exx(macroElementProps.yPos,macroElementProps.xPos)=Exx;
        DV.Eyy(macroElementProps.yPos,macroElementProps.xPos)=Eyy;
        DV.t(macroElementProps.yPos,macroElementProps.xPos)=Theta;
        DV.w(macroElementProps.yPos,macroElementProps.xPos)=v;
        
        %                 objectiveValue = ObjectiveCalculateEffectiveVars(x,DmatrixIN, matProp,config);
        
        rhoArray = [rhoArray;v];
        ExxArray=[ExxArray;Exx];
        EyyArray=[EyyArray;Eyy];
        thetaArray=[thetaArray;Theta];
    else
        DV.Exx(macroElementProps.yPos,macroElementProps.xPos)=ActualExx;
        DV.Eyy(macroElementProps.yPos,macroElementProps.xPos)=ActualEyy;
        DV.t(macroElementProps.yPos,macroElementProps.xPos)=ActualThetaValue;
        DV.w(macroElementProps.yPos,macroElementProps.xPos)=0;
        
    end
end

% save the Exx field
outname = sprintf('./out%i/ExxSubSysValues%i.csv',folderNum,macro_meso_iteration);
csvwrite( outname,DV.Exx);

% save the Eyy field
outname = sprintf('./out%i/EyySubSysValues%i.csv',folderNum,macro_meso_iteration);
csvwrite( outname,DV.Eyy);

% save the Theta field
outname = sprintf('./out%i/ThetaSubSysValues%i.csv',folderNum,macro_meso_iteration);
csvwrite( outname, DV.t);

% save the density field
outname = sprintf('./out%i/densityUsedSubSysValues%i.csv',folderNum,macro_meso_iteration);
csvwrite( outname,  DV.w);

figure(1)
if(1==2)
    scatter3(ExxArray,EyyArray,rhoArray);
    xlabel('Exx');
    ylabel('Eyy');
    zlabel('Rho,Density');
    
    figure(2)
    sf = fit([ExxArray, EyyArray],rhoArray,'poly23')
    plot(sf,[ExxArray,EyyArray],rhoArray)
    
    xlabel('Exx');
    ylabel('Eyy');
    zlabel('Rho,Density');
    
    figure(3)
    scatter3(ExxArray,EyyArray,thetaArray);
    xlabel('Exx');
    ylabel('Eyy');
    zlabel('theta');
    
    figure(4)
    histogram(thetaArray)
    
%     figure(5)
%     DV = DV.CalculateVolumeFractions( config,matProp);
else
    DV.CalculateVolumeFractions( config,matProp)
    p = plotResults;
    FEACalls=1;
    p.plotTopAndFraction(DV,  config, matProp, FEACalls); % plot the results.
    
   
    nameGraph = sprintf('./MesoDesignExxEyyThetaVars%i.png', config.macro_meso_iteration);
    print(nameGraph,'-dpng');
end
% 
 p = plotResults;
diffExx = ExxMacro- DV.Exx;
diffEyy = EyyMacro- DV.Eyy;
diffTheta = ThetaMacro- DV.t;


relativeErrorExx=diffExx./ExxMacro;
relativeErrorEyy=diffEyy./EyyMacro;
relativeErrorTheta=diffTheta./ThetaMacro;


% make the range -1 to 1
relativeErrorExx(relativeErrorExx>1)=1;
relativeErrorExx(relativeErrorExx<-1)=-1;

relativeErrorEyy(relativeErrorEyy>1)=1;
relativeErrorEyy(relativeErrorEyy<-1)=-1;

relativeErrorTheta(relativeErrorTheta>1)=1;
relativeErrorTheta(relativeErrorTheta<-1)=-1;

% figure
% p.PlotArrayGeneric( diffExx, 'diffExx')
% figure
% p.PlotArrayGeneric( diffEyy, 'diffEyy')
% figure
% p.PlotArrayGeneric( diffTheta, 'diffTheta')

figure
subplot(2,2,1)
p.PlotArrayGeneric(100* relativeErrorExx, 'Percent Error Exx')
subplot(2,2,2)
p.PlotArrayGeneric( 100*relativeErrorEyy, 'Perecent Error Eyy')
subplot(2,2,3)
p.PlotArrayGeneric(100* relativeErrorTheta, 'Percent Error Theta')
subplot(2,2,4)
p.PlotArrayGeneric(diffTheta, 'Diff Theta')
nameGraph = sprintf('./MesoDesignExxEyyThetaVarsPercentError%i.png', config.macro_meso_iteration);
print(nameGraph,'-dpng');



