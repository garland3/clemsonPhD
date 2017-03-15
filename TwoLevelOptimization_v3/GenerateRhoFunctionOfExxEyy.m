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
%                   outname = sprintf('./out%i/DsystemIter%i_Element_%i.csv',folderNum,macro_meso_iteration,elementNumber);
        Din = csvread(outname);
        
        
        mode = 1;
        doPlot=0;
        
        
        minValueForZeroTerms = 1000000;
        thetaOfMin = 0;
        count = 1;
        
        arrayOfOuts = [];
        %             thetaValuesToTest=-pi/2:pi/360:pi/2;
        temp33=pi/2-ActualThetaValue;
        incr = pi/180;
        thetaValuesToTest1=temp33-4*pi/8+incr:incr:temp33+4*pi/8-incr;
        thetaValuesToTest=-ActualThetaValue+incr:incr:pi-ActualThetaValue-incr;
        %                 temp0=100;
        
        
        temp0=100;
        closestTheta = 0;
        
        for theta = thetaValuesToTest
            Dout = matProp.rotateDmatrix(config,theta, Din);
            sumZeroTerms =abs( Dout(1,3))+abs(Dout(2,3))+abs(Dout(3,1))+abs(Dout(3,2));
            
            
            
            if(sumZeroTerms<minValueForZeroTerms)
                thetaOfMin=theta;
                minValueForZeroTerms=sumZeroTerms;
            end
            arrayOfOuts=[arrayOfOuts;sumZeroTerms];
            
            diffTheta = abs(theta-ActualThetaValue);
            if(diffTheta<temp0)
                closestTheta=theta;
                temp0=diffTheta;
                
            end
            
            count=count+1;
        end
        Theta=thetaOfMin;
        
         if(Theta<0)
            Theta=Theta+pi/2;
        end
        if(Theta>pi/2)
            Theta=Theta-pi/2;
        end
        
        if doPlot==1
            
            %              actualTheta = zeros(size(thetaValuesToTest));
            %                 xThetasActual = [thetaValuesToTest(1) pi/2-closestTheta thetaValuesToTest(end)];
            %                 xpropsedTheta = [thetaValuesToTest(1) Theta thetaValuesToTest(end)];
            xThetasActual = [thetaValuesToTest(1) closestTheta thetaValuesToTest(end)];
            xpropsedTheta = [thetaValuesToTest(1) pi/2-Theta thetaValuesToTest(end)];
            yThetas = [0 1 1]*matProp.E_material1;
            
            plot(thetaValuesToTest,arrayOfOuts)
            hold on
            stairs(xThetasActual,yThetas,'g-o');
            stairs(xpropsedTheta,yThetas,'r-x');
            title(sprintf('e=%i,Minimiz terms, Green = actual theta, Red = proposed Theta',e));
            
            hold off
        end
        
        
       
        
        %             cases = [Theta pi/2-Theta ];
        %             maxEs=0;
        denominator = 1-matProp.v;
        %              ExxFinal =0;
        %              EyyFinal=0;
        %              thetaFinal=0;
        %             for mmm = cases
        Dout = matProp.rotateDmatrix(config,Theta, Din);
        Exx=Dout(1,1)*denominator;
        Eyy=Dout(2,2)*denominator;
        Theta = pi/2-Theta;
        
        % Scale UP based on the rho (x) density
        Exx=Exx*1/(macroElementProps.densitySIMP);
        Eyy=Eyy*1/(macroElementProps.densitySIMP);
        
        if(ActualExx>ActualEyy)
            if(Exx>Eyy)
                % great, do nothing
            else
                % otherwise, flip
                temp= Exx;
                Exx=Eyy;
                Eyy=temp;
            end
        else %ActualExx<ActualEyy
            if(Exx<Eyy)
                % great, do nothing
            else
                % otherwise, flip
                temp= Exx;
                Exx=Eyy;
                Eyy=temp;
            end
            
        end
        
        diffX = ActualExx-Exx;
        diffY = ActualEyy-Eyy;
        totalDiff = diffX+diffY;
        
        
        
        
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
        DV.Exx(macroElementProps.yPos,macroElementProps.xPos)=matProp.E_material2/2;
        DV.Eyy(macroElementProps.yPos,macroElementProps.xPos)=matProp.E_material2/2;
        DV.t(macroElementProps.yPos,macroElementProps.xPos)=0;
        DV.w(macroElementProps.yPos,macroElementProps.xPos)=0;
        
    end
end

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
    
    figure(5)
    DV = DV.CalculateVolumeFractions( config,matProp);
else
    DV.CalculateVolumeFractions( config,matProp)
    p = plotResults;
    FEACalls=1;
    p.plotTopAndFraction(DV,  config, matProp, FEACalls); % plot the results.
    
    loadcaseIndex=1;
    nameGraph = sprintf('./MesoDesignExxEyyThetaVars%i.png', config.macro_meso_iteration);
    print(nameGraph,'-dpng');
end
