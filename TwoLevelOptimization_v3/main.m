function main()
clear
% 
%     combinedTopologyOptimization('1', '1', '1','60', 'na');
%     combinedTopologyOptimization('1', '1', '1','200', 'na');

    
 combinedTopologyOptimization('1', '1', '1','202', 'na');

 
% test meso design for element 1, macro iteration 1
%  combinedTopologyOptimization('1', '1', '1','100', '781');
% combinedTopologyOptimization('1', '1', '1','100', '191');
%   combinedTopologyOptimization('1', '1', '1','100', '2');
%     combinedTopologyOptimization('1', '1', '1','100', '252');
%       combinedTopologyOptimization('1', '1', '1','100', '721');
 
 

% numloops = 5;
% for i =1:numloops
% %     combinedTopologyOptimization('1', 'not used', 1,i)
% %     combinedTopologyOptimization(useInputArgs, w1text, macro_meso_iteration,mode, singleMeso_elementNumber);
%      combinedTopologyOptimization('no', 'NA', i,'na', 'na');
% end
% close all
% 
% 
%   combinedTopologyOptimization('1', '1',numloops,'12', 'na');

% folderNum = 0;
% all = [];
% iterationDiff = [];
% for i =1:numloops
%     macro_meso_iteration = i;
%     outname = sprintf('./out%i/storeOptimizationVarMacroLoop%i.csv',folderNum,macro_meso_iteration);
%     storeOptimizationVar=  csvread(outname);
%     all=[all;storeOptimizationVar];
%     [rr,cc] = size (storeOptimizationVar);
%     iterationDiff= [iterationDiff zeros(1,rr-1)];
%     iterationDiff = [iterationDiff 1];
% end
% designVars = DesignVars;
% designVars.storeOptimizationVar=all;
% p = plotResults;
% [r,c] = size (all);
% loopNumb=r;
% 
% % Get size of the grid
% previousIterationNum=1;
% outname = sprintf('./out%i/volfractionfield%i.csv',folderNum,previousIterationNum);
% [gridrows, gridcolumns]=size(csvread(outname));
% totalvolume = gridrows*gridcolumns;
% 
% % normalize the volumes
% %mtemp = max(designVars.storeOptimizationVar(1:loopNumb,4));
% % designVars.storeOptimizationVar(1:loopNumb,4) = designVars.storeOptimizationVar(1:loopNumb,4)/totalvolume;
% % %mtemp2 = max(designVars.storeOptimizationVar(1:loopNumb,5));
% % designVars.storeOptimizationVar(1:loopNumb,5) = designVars.storeOptimizationVar(1:loopNumb,5)/totalvolume;
% 
% 
% 
% settings = 0;
% matProp = 0;
% hold on
% p.PlotDesignMetrics( designVars, settings, matProp, loopNumb)
% xvalues = 1:loopNumb;
% stairs(xvalues,iterationDiff);
% hold off
% 
% w1 = 0;
%  nameGraph = sprintf('./optiParaViaIterations%i.png', w1);
% print(nameGraph,'-dpng');


end