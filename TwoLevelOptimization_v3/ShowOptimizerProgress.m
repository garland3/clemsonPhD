% --------------------------------------------------------
% --------------------------------------------------------
% -----------------------
% SHOW STATUS OF OPTIMIZER
% Print and Plot information.
% Store data in the storeOptimizationVar
% -----------------------
% --------------------------------------------------------
% --------------------------------------------------------
function ShowOptimizerProgress(DV,doPlot,name,FEACalls,config, matProp)


% PRINT RESULTS
disp([' FEA calls.: ' sprintf('%4i',FEACalls) ' Obj.: ' sprintf('%10.4f',DV.c) ' Vol. 1: ' sprintf('%6.3f', DV.currentVol1Fraction)  ' Vol. 2: ' sprintf('%6.3f', DV.currentVol2Fraction) ...
    ' Target E.: ' sprintf('%4i',config.targetAvgExxEyy)    ' Current E.: ' sprintf('%4i',DV.actualAverageE) ' Avg Meso Density '  sprintf('%4i',DV.averageMesoDensity) name] );



if(doPlot==1 &&  config.RunPalmetto==0)
    p = plotResults;
    p.plotTopAndFraction(DV,  config, matProp, FEACalls); % plot the results.
end
