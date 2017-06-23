% --------------------------------------------------------
% --------------------------------------------------------
%
% test for termaination of the function.
% normalize the objectives, compare to a moving average. If moving average
% below target, then return 1
% if not terminated, then return 0
% --------------------------------------------------------
% --------------------------------------------------------
function status = TestForTermaination(DV, config,masterloop)

if(masterloop<3)
    status=0;
    return;
end
status = 0;
y2 = DV.storeOptimizationVar(:,2); % Elastic Compliance

t = config.terminationAverageCount;
if(size(y2)<(t+3))
    return;
end


y3 = DV.storeOptimizationVar(:,3); % Heat Compliance
y4 = DV.storeOptimizationVar(:,4); % material 1
y5 = DV.storeOptimizationVar(:,5); % material 2

avg_y2 = FindAvergeChangeOfLastValue(y2,config); % elastic objective
avg_y3 = FindAvergeChangeOfLastValue(y3,config); % heat objective
avg_y4 = FindAvergeChangeOfLastValue(y4,config); % material 1
avg_y5 = FindAvergeChangeOfLastValue(y5,config); % material 2




tt = config.terminationCriteria;
if(        avg_y2<tt ...
        && avg_y3<tt ...
        && avg_y4<tt ...
        && avg_y5<tt)
    status=1;
    % print the vars to screen
    [avg_y2 avg_y3 avg_y4 avg_y5 ]
end