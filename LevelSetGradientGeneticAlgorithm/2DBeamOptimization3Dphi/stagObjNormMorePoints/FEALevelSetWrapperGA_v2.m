function [obj, cinequality] = FEALevelSetWrapperGA_v2(zpoints)

[xpoints,ypoints] = meshgrid(0:0.25:6,0:0.25:2); 
zpoints_v2 = xpoints;
[xdim, ydim] = size(zpoints_v2);

if(zpoints ==1.111)
    zpoints = ones(1,xdim*ydim)*1.1;
     zpoints = [0.77016,0.79896,0.9748,0.96156,1,1,1,1,0.98796,1,1,0.81731,1,0.75,0.75825,1,0.71931,0.75,0.31351,0.75,1,0.21805,0.69285,0.3325,0.86142,0.5,0.5625,0.97387,1,0.77509,1,0.0625,0.125,0.39948,0.25,0.5,0.3125,1,0.13457,0.46967,0.1569,0.19179,0.0625,1,0.34246,0.091252,1,0.17655,0.9257,0,0.10174,0,0.0625,1,0.20823,0.038071,1,0.21263,0.097803,0.02318,0.75,0.13212,0.17747,1,0.11177,0.25,0.06737,0.83458,0,0.73948,0.625,0.11428,1,0.27905,0.0625,0.25,0.060259,0.0533,0.8066,0.056521,0.20952,0.14371,0.88298,0.13154,1,0.011208,0.16171,0.12869,0.71867,0.074246,0.48458,0.44058,0.40159,0.56905,0.11839,0.89883,0.52232,0.75,0.099992,0.22677,0.21971,0.1605,1,0.16744,0.16426,0.16338,0.030085,0.32954,1,0.18605,0.1356,0.053788,0.18368,1,0.45246,0.31865,0.35013,0.12168,0.76253,1,0.55992,0.62733,0.15619,0.34872,0.29318,0.45233,0.55019,0.49973,0.10064,0.16351,0.52764,0.17871,0.47611,0.19235,0.91321,0.22533,0.47556,0.102,0.36616,1,0.1192,0.22033,1,0.52446,0.1875,0.099571,0.16465,0.32272,0.3288,0.38825,0.027709,0.67219,0.027364,0,0.0625,0.23719,0.068477,0.25,0.14384,0.81139,0.1021,0.20523,0.34451,0.69744,1,1,1,0.024471,0.066663,0.40403,0.25,0.21695,0.20792,0.13632,0.34491,0.89066,0.81374,1,1,0.93767,1,1,0.74166,0.94633,1,0.99807,0.7107,0.75,0.4251,0.27311,0.28021,0.027986,0.25,0.34954,0.36167,0.045209,0.30321,0.29987,0.26748,0.0625,0.23772,0.5,0.24122,0.19578,0.5,0.16153,0.17801,0,0.455,0.25,0.0625,0.0625,0,0.43313,0.13025,0.52349,0.14303,0.23892,0.24526,0.20834,0.25152,0.23448,0.3657,0.086424,0.27507]
    %zpoints = [0.80419,0.92926,0.9231,0.92713,0.82608,0.82837,0.96309,0.27148,0.43926,0.99629,0.24976,0.30176,0.26478,0.05422,0.89502,0.76891,0.44518,0.047689,0.32923,0.062045,0.10381]

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

weightCost = 0.7;
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




