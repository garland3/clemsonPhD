function [obj, cinequality] = FEALevelSetWrapperGA_v2(zpoints)

[xpoints,ypoints] = meshgrid(0:0.5:6,0:0.5:2); 
zpoints_v2 = xpoints;
[xdim, ydim] = size(zpoints_v2);


if(zpoints ==1)
    zpoints = ones(1,xdim*ydim)*1.1;
    % zpoints = [0.53767,1,-1,0.86217,0.31877,-1,-0.43359,0.34262,1,1,-1,1,0.7254,-0.063055,0.71474,-0.20497,-0.12414,1,1,1,0.6715,-1,0.71724,1,0.48889,1,0.72689,-0.30344,0.29387,-0.78728,0.8884,-1,-1,-0.8095,-1,1,0.32519,-0.75493,1,-1,-0.10224,-0.24145,0.31921,0.31286,-0.86488,-0.030051,-0.16488,0.62771,1,1,-0.86365,0.077359,-1,-1,-0.0068493,1,-0.76967,0.37138,-0.22558,1,-1,0.032557,0.55253,1,1]
    zpoints = [ 0.59707,0.57584,0.1,0.12539,0.1875,0.56948,0.45,0.1,0.45,0.7,0.45,0.70327,-0.05,-0.4,-0.50036,0.1,0.1,0.11563,-0.82878,-0.19513,0.18995,-0.97804,-0.89227,0.32353,-0.2,-0.001341,0.52886,-0.089062,0.1,-0.5563,-0.9,-0.62044,-0.8141,0.1,0.6,0.6,-0.0053513,-0.99429,-0.4,0.77578,0.93438,0.17813,-0.15,0.16823,0.39297,0.6,0.95,-0.64902,-0.82018,-0.39414,-0.45223,-0.93101,-0.05,-0.55,0.60317,0.95,-0.040442,-0.69433,-0.85604,-0.9,-0.38438,-0.49533,-0.81483,-0.9,-0.4]
    doplot = 1;
    
    
else
    doplot = 0;
end


count = 1;
for i = 1:xdim
    zpoints_v2(i,:) = zpoints(count:ydim+count-1);
    count = count + ydim;  
end




[maxVonMisesStress, cost,maxDisplacement] = FEALevelSet_2D_v8(xpoints,ypoints,zpoints_v2, doplot); % Call FEA for first time

% vonMises Stress
% 585.1768 = solid = min, 
% 3000=  max
normalizedStress = (maxVonMisesStress-585)/(3000-585);
normalizedCost = cost;

weightCost = 0.5;
weightStress = 0.5;

obj = normalizedStress*weightStress  +   normalizedCost*weightCost % mininimize vonMisesStress, 
cinequality = [normalizedStress-1; maxDisplacement-0.3] ;  % sstress must be below 3000, max dipslacement below 0.3
% cequality = ydisplacment+0.0623; % displacement will be negative, so add 0.0623 to make = 0
%feature('GetPid')





