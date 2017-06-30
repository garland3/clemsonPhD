close all
MacroExxColumnTotal=[];
MacroEyyColumnTotal=[];
MacroThetaColumnTotal=[];
ErrorData=[];

% combinedTopologyOptimization('1', '1', '5','203', 'na');

folderName='data';

useSubSysValues = 0;

macro_meso_iteration=5;
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


size(MacroExxColumnTotal)

config = Configuration;
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
numBins = 30;
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






f1 = fit([x y],z,'poly22',options)
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
%     zlim([0 1])
%     size(RhoColumnTotal)