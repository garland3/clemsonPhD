function [] = plotMode90Data(iterations)

macroModel=[];
actualModel=[];
density=[];
for macro_meso_iteration = 1:iterations
    outname = sprintf('./mode90/mode90_objective%i.csv',macro_meso_iteration);
    temp=csvread(outname);
    actualModel=[actualModel temp];
    
    
    outname = sprintf('./out%i/storeOptimizationVarMacroLoop%i.csv',0,macro_meso_iteration);
    storeOptimizationVar=csvread(outname);
    elasticObjectiveMacroModel = storeOptimizationVar(end,2); % Elastic Compliance
    macroModel=[macroModel elasticObjectiveMacroModel];
    
    
    outname = sprintf('./mode90/Density%i.csv',macro_meso_iteration);
    densityTemp=csvread(outname);
    density=[density densityTemp];
    
end

actualModel
macroModel
density
x=1:size(actualModel,2)
plot(x,actualModel,'b-*')
hold on
plot(x,macroModel,'r-*')
plot(x,density,'g-*');
title('Macro Model elastic value in red, Elastic value based on using D_h in blue, Density in Green')
hold off
% fprintf('Mode 90, final objective value uisng D_sub values %f, compared to %f using macro model\n',DV.c,elasticObjectiveMacroModel);
%     plotMode90Data(config.macro_meso_iteration);

nameGraph = sprintf('./mode90_iteration%i.png',macro_meso_iteration);
print(nameGraph,'-dpng');

end