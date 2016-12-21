clear
clc
close all

% material properties Object
matProp = MaterialProperties;
settings = Configuration;
plotter = plotResults;
designVars = DesignVars(settings);

endingstoreOptimizationVar = []

w1 = 0.5;
if(w1 ==1)
    outname11 = sprintf('./temp/out11/storeOptimizationVar.csv');
    outname12 = sprintf('./temp/out12/storeOptimizationVar.csv');
    outname13 = sprintf('./temp/out13/storeOptimizationVar.csv');
    outnameGradMat =  sprintf('./out10/storeOptimizationVar.csv');
elseif(w1 == 0.5)
%     outname11 = sprintf('./temp/out14/storeOptimizationVar.csv');
%     outname12 = sprintf('./temp/out15/storeOptimizationVar.csv');
%     outname13 = sprintf('./temp/out16/storeOptimizationVar.csv');
    outname11 = sprintf('./out14/storeOptimizationVar.csv');
    outname12 = sprintf('./out15/storeOptimizationVar.csv');
    outname13 = sprintf('./out16/storeOptimizationVar.csv');
    outnameGradMat =  sprintf('./out5/storeOptimizationVar.csv');
end

filenames = {outname11 outname12 outname13 outnameGradMat}
for i = 1:4
    fname = char(filenames(i));
    designVars.storeOptimizationVar = csvread(fname);
    
    
     % get information about the final designs and store it,
    last = size(designVars.storeOptimizationVar,1);
    temp = [designVars.storeOptimizationVar(last,:)];
    endingstoreOptimizationVar = [endingstoreOptimizationVar; temp ];
end

% designVars.storeOptimizationVar = [designVars.storeOptimizationVar;designVars.c, designVars.cCompliance, designVars.cHeat,vol1Fraction,vol2Fraction,fractionCurrent_V1Local,densitySum];

complianceObj = endingstoreOptimizationVar(:,2)
heatObj = endingstoreOptimizationVar(:,3)
vol1Const = endingstoreOptimizationVar(:,4);
vol2Const = endingstoreOptimizationVar(:,5);

complianceObj = complianceObj/max(complianceObj);
heatObj = heatObj/max(heatObj);



% 

if(w1 ==1)
    yBar = [complianceObj  vol1Const vol2Const];
elseif(w1 ==0.5)
    yBar = [complianceObj heatObj vol1Const vol2Const];
end

b = bar(yBar);

set(gca,'XTickLabel',{'single Mat 1', 'single Mat 2', 'single Mat 0.5', 'Gradient'})

if(w1 ==1)
    legend(b,{'elastic','vol 1', 'vol 2'})
    print('matComparisonW1.png','-dpng')
elseif(w1 ==0.5)
    legend(b,{'elastic','heat','vol 1', 'vol 2'})
     print('matComparisonW0.5.png','-dpng')
end

csvwrite('yBarFromPlotting.csv',yBar)