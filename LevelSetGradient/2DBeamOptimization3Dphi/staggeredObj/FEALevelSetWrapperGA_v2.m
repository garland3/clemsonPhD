function [obj, cinequality] = FEALevelSetWrapperGA_v2(zpoints)

[xpoints,ypoints] = meshgrid(0:0.5:6,0:0.5:2); 
zpoints_v2 = xpoints;
[xdim, ydim] = size(zpoints_v2);

if(zpoints ==1.111)
    zpoints = ones(1,xdim*ydim)*1.1;
    zpoints = [0.41326,0.91945,0.32528,1.0766,0.7604,0.59246,0.53448,-0.12955,-1.5242,-0.56249,-0.53109,-1.9071,1.3655,-1.0154,-1.7826,-1.681,-1.8103,-1.9088,-1.3607,-1.4429,1.9716,1.9752,0.82417,0.42362,0.50392,-0.31075,-0.82959,-1.8288,-1.5034,-1.3352,0.43713,0.014995,0.28565,-0.71902,-1.9437,-0.73974,-0.86812,-1.5408,-0.76833,0.33211,0.46958,0.06198,0.50054,-0.45375,-0.57695,-1.9792,-1.7148,-0.64697,-1.9214,-0.30002,-0.4617,-0.23326,-0.11518,-1.4677,-1.1271,-1.8548,-1.3892,-1.4264,-0.88006,-0.72142,-0.94898,-0.53901,-1.2024,-0.020559,-0.80791]

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

weightCost = 0.5;
weightStrainE = 1-weightCost;


obj = (strainEnergyNormalized*weightStrainE ...
    + weightCost*cost)*1000; % do a weighted sum. Multiply by 1000 to scale up so we don't get close the the ending conditions for the GA. 

targetCost = 0.75;

% if the cost goes below, 0.5, then just leave the objective calc like it
% is at 0.5. 
if(cost<=targetCost)
    obj = (strainEnergyNormalized*weightStrainE ...
    + weightCost*targetCost)*1000;
end

cinequality = [cost] ;  % sstress must be below 3000, max dipslacement below 0.3


fprintf('Strain E  = %f and cost = %f  and obj = %f\n',strainEnergy, cost, obj)




