function [obj, cinequality] = FEALevelSetWrapperGA_v2(zpoints)

[xpoints,ypoints] = meshgrid(0:0.5:6,0:0.5:2); 
zpoints_v2 = xpoints;
[xdim, ydim] = size(zpoints_v2);

if(zpoints ==1.111)
    zpoints = ones(1,xdim*ydim)*1.1;
    zpoints = [0.18912,0.22217,0.29169,0.96246,0.89994,0.29282,0.25598,0.22466,0.21207,0.10863,0.99593,0.64832,0.68755,0.98447,0.92239,0.85693,0.042297,0.15389,0.94481,0.96876,0.96751,0.96297,0.92239,0.57985,0.15883,0.10416,0.049583,0.061613,0.052708,0.23621,0.71678,0.61582,0.42264,0.29713,0.11235,0.076414,0.090168,0.30447,0.47919,0.29631,0.4619,0.52948,0.56037,0.41421,0.16352,0.33169,0.17031,0.47967,0.39828,0.45968,0.36776,0.18924,0.80956,0.50095,0.24261,0.055334,0.35118,0.22035,0.42321,0.15058,0.01708,0.38557,0.22433,0.41518,0.29625]

    doplot = 1;
else
    doplot = 0;
end

zpoints = zpoints*4-2; % un-normalize the values. 

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




