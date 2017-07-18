close all
MacroExxColumnTotal=[];
MacroEyyColumnTotal=[];
MacroThetaColumnTotal=[];
ErrorData=[];

% combinedTopologyOptimization('1', '1', '5','203', 'na');

folderName='data';

useSubSysValues = 0;

macro_meso_iteration=1;
folderName='ErrorData';
outnameExx = sprintf('./%s/MacroExxColumnTotal%i.csv',folderName,macro_meso_iteration);

% outnameEyy = sprintf('./%s/MacroEyyColumn%i.csv',folderName,macro_meso_iteration);
outnameEyy = sprintf('./%s/MacroEyyColumnTotal%i.csv',folderName,macro_meso_iteration);


% outnametheta = sprintf('./%s/MacroThetaColumn%i.csv',folderName,macro_meso_iteration);
outnametheta = sprintf('./%s/MacroThetaColumnTotal%i.csv',folderName,macro_meso_iteration);

outnameError = sprintf('./%s/ErrorTotal%i.csv',folderName,macro_meso_iteration);
% outnamerho = sprintf('./%s/RhoColumn%i.csv',folderName,macro_meso_iteration);




%------------------------
% read the macro columns as well. This will help with future analysis
% ----------------------------
% save the MacroExxColumn
MacroExxColumn=csvread(outnameExx);
MacroExxColumnTotal=[MacroExxColumnTotal; MacroExxColumn];

% save the MacroEyyColumn
temp=csvread(outnameEyy);
MacroEyyColumnTotal=[MacroEyyColumnTotal; temp];

% save the MacroThetaColumn%
temp=csvread(outnametheta);
MacroThetaColumnTotal=[MacroThetaColumnTotal; temp];

% save the Error
temp=csvread(outnameError);
ErrorData=[ErrorData; temp];


% Select a subset
step=11;
temp=config.validationGridSizeNelx-1;
thetaVector = 0:(pi/2)/temp:pi/2;


Down = thetaVector(step)-0.001
Up=thetaVector(step)+0.001
logicTest1=Down<MacroThetaColumnTotal  ;
logicTest2=MacroThetaColumnTotal<Up;
logicTest=(logicTest1+logicTest2)>1.1

MacroExxColumnTotal=MacroExxColumnTotal(logicTest);
MacroEyyColumnTotal=MacroEyyColumnTotal(logicTest);
MacroThetaColumnTotal=MacroThetaColumnTotal(logicTest)
ErrorData=ErrorData(logicTest);


size(MacroExxColumnTotal)

config = Configuration;

if(1==1)
    % if(config.UseLookUpTableForPsuedoStrain==1)
    %     m=config.lookupSearchScheme;
    % else
    %
    %      m = config.mesoVolumeUpdateMethod;
    % end
    % m=2
    % outname = sprintf('./%s/Valid_TotalError_Method%i_Config%i.csv',folderName,config.UseLookUpTableForPsuedoStrain,m);
    % totalErrorV2 =     csvread(outname);
    
    options= fitoptions;
    % %    options.Normalize ='on';
    % %    options.fittype='poly22';
    
    
    
    [numDataPoints]=size(ErrorData,1);
    
    % ---------------------------------
    %
    %
    %            Average the data into bins
    %
    %
    % ---------------------------------
    numBins = 50;
    % [xbins, ybins, zbins]=meshgrid(1:numBins,1:numBins,1:numBins);
    [xbins, ybins]=meshgrid(1:numBins,1:numBins);
    totalBins = numBins^2 % only consider the Exx, And Eyy
    SumInBin=xbins*0;
    CountInBin=SumInBin;
    
    matProp = MaterialProperties; % material properties Object
    step = matProp.E_material1/(numBins);
    % stepTheta = (pi/2)/numBins;
    for lll = 1:totalBins
        xIndex = xbins(lll);
        xLowerLimit = step*(xIndex-1);
        xUpperLimit = step*(xIndex);
        
        yIndex = ybins(lll);
        yLowerLimit = step*(yIndex-1);
        yUpperLimit = step*(yIndex);
        
        
        %     zIndex = zbins(lll);
        %     xLowerLimit = stepTheta*(xIndex-1);
        %     xUpperLimit = stepTheta*(xIndex);
        
        
        for kk = 1:numDataPoints
            ExxValue = MacroExxColumnTotal(kk);
            EyyValue = MacroEyyColumnTotal(kk);
            %           thetaValue = MacroExxColumnTotal(kk);
            
            if(xLowerLimit<ExxValue && ExxValue<=xUpperLimit)
                if(yLowerLimit<EyyValue && EyyValue<=yUpperLimit)
                    errorValue = ErrorData(kk);
                    
                    errorValue=min(errorValue,3);
                    
                    
                    SumInBin(lll)= SumInBin(lll)+errorValue;
                    CountInBin(lll)= CountInBin(lll)+1;
                end
            end
        end
    end
    
    %
    % x =MacroExxColumnTotal;
    % y = MacroEyyColumnTotal;
    % z=ErrorData;
    
    newExxArray=[];
    newEyyArray=[];
    newErrorArray=[];
    
    for lll = 1:totalBins
        c= CountInBin(lll);
        if(c>0)
            xIndex = xbins(lll);
            xLowerLimit = step*(xIndex-1);
            xUpperLimit = step*(xIndex);
            xavg=(xLowerLimit+xUpperLimit)/2;
            
            yIndex = ybins(lll);
            yLowerLimit = step*(yIndex-1);
            yUpperLimit = step*(yIndex);
            yavg=(yLowerLimit+yUpperLimit)/2;
            
            avgError = SumInBin(lll)/CountInBin(lll);
            
            newExxArray=[newExxArray;xavg];
            newEyyArray=[newEyyArray;yavg];
            
            newErrorArray=[newErrorArray ;avgError];
        end
    end
    
    x =newExxArray;
    y = newEyyArray;
    z=newErrorArray;
    
    
    
    
    
    
    f1 = fit([x y],z,'poly33',options)
    coeffvalues(f1)
    
    figure
    plot(f1, [x y], z);
    hold on
    scatter3(x,y,z,'b')
    % %     hold on
    % %     scatter3(x,y,annZ,'r');
    
    title('Fit with data points.  Blue=Actual Error ')
    xlabel('Exx');
    ylabel('Eyy');
    zlabel('error');
    
elseif(1==0)
    
    setdemorandstream(491218382)
    % net = fitnet(10,'trainbfg');
    %  net = cascadeforwardnet(10);
    net = feedforwardnet(10);
    figure(1)
    % view(net)
    % print('ANN_view.png','-dpng');
    
%     MacroExxColumnTotal=[MacroExxColumnTotal; MacroExxColumn];
% 
% % save the MacroEyyColumn
% temp=csvread(outnameEyy);
% MacroEyyColumnTotal=[MacroEyyColumnTotal; temp];
% 
% % save the MacroThetaColumn%
% temp=csvread(outnametheta);
% MacroThetaColumnTotal=[MacroThetaColumnTotal; temp];
% 
% % save the Error
% temp=csvread(outnameError);
% ErrorData=[ErrorData; temp];
    
    x=[MacroExxColumnTotal';  MacroEyyColumnTotal';  MacroThetaColumnTotal'];
    t = [ErrorData'];
    
    size(x)
    size(t)
    
    % Train
    [net,tr] = train(net,x,t);
    
    plotperform(tr)
    print('ANN_preformance.png','-dpng');
    
    
    y = net(x);
    plotregression(t,y)
    print('ANN_regressionTest2.png','-dpng');
    
end
%     zlim([0 1])
%     size(RhoColumnTotal)