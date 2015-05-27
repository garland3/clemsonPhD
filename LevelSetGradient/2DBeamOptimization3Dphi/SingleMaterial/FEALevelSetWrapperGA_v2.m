function [obj, cinequality] = FEALevelSetWrapperGA_v2(zpoints)

[xpoints,ypoints] = meshgrid(0:0.5:6,0:0.5:2); 
zpoints_v2 = xpoints;
[xdim, ydim] = size(zpoints_v2);

if(zpoints ==1.111)
    zpoints = ones(1,xdim*ydim)*1.1;
    zpoints = [1.1168,0.56374,1.0158,0.88715,0.98241,0.52183,-1.2413e-05,-0.60339,1.4242,0.45132,-1.9005,1.3065,-0.85271,-1.0997,-1.0787,-1.0484,-1.1419,-1.0327,-1.1223,1.3372,1.5049,-1.8732,-0.17637,1.1222,-0.52332,-1.5837,-0.13839,-0.35297,-1.5033,-0.69575,1.2737,1.7734,-1.3552,-1.8069,-0.46444,0.25,-0.78003,-2,-0.39282,-0.47309,-0.88808,0.76532,0.35769,-0.68493,-0.37962,0.76413,0.00011905,0.88513,-1.8164,-0.32741,-1.9826,-0.24331,0.28242,1.3841,-1.0661,-1.4999,-1.5946,-1.4629,-1.443,-1.9713,-1.7268,-1.0931,-0.8182,-1.8273,-0.87771]
    doplot = 1;
else
    doplot = 0;
end


count = 1;
for i = 1:xdim
    zpoints_v2(i,:) = zpoints(count:ydim+count-1);
    count = count + ydim;  
end


[maxVonMisesStress, cost,maxDisplacement,strainEnergy] = FEALevelSet_2D_v8(xpoints,ypoints,zpoints_v2, doplot); % Call FEA for first time

minStrainE = 40;
maxStrainE = 500;
strainEnergyNormalized = (strainEnergy-minStrainE)/(maxStrainE-minStrainE);

weightCost = 0.9;
weightStrainE = 1-weightCost;


obj = (strainEnergyNormalized*weightStrainE ...
    + weightCost*cost)*1000; % do a weighted sum. Multiply by 1000 to scale up so we don't get close the the ending conditions for the GA. 

targetCost = 0.3;

% if the cost goes below, 0.5, then just leave the objective calc like it
% is at 0.5. 
if(cost<=targetCost)
    obj = (strainEnergyNormalized*weightStrainE ...
    + weightCost*targetCost)*1000;
end

cinequality = [cost] ;  % sstress must be below 3000, max dipslacement below 0.3


fprintf('Strain E  = %f and cost = %f  and obj = %f\n',strainEnergy, cost, obj)




