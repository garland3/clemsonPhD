function [obj, cinequality, cequality] = FEALevelSetWrapperGA(zpoints)



[xpoints,ypoints] = meshgrid(0:1:10,0:0.25:1);
zpoints_v2 = zeros(5,11);

% Wrapper version 5
% FEA version 7

% after 2500 evaluations and 50 generations
% cost = 0.5173
%zpoints = [0.10689	0.64314	0.70731	-0.27855	-1.5	-1.31067	-0.45423	-0.52377	-0.43743	-1.03323	-0.2049	-0.04854	-0.01329	-0.18558	0.55965	0.62553	0.18225	-0.31614	-0.30015	-1.02867	-0.53433	-0.88878	-1.49964	1.04412	0.78072	-0.38154	0.04917	-0.33039	-0.02682	-0.33582	0.57312	-1.38972	0.60729	0.18489	1.27056	0.4704	0.31509	0.47898	-0.26517	-0.85263	-1.44156	-1.48608	-1.27854	-0.52986	0.23031	-1.27836	-0.05139	-1.43739	-1.39227	0.82536	-0.83418	0.94839	-0.98382	-1.15284	0.54333]

% using the best 45 of the previous optimization as the starting points, I
% continued the optimization with another 2500 evalutions
% zpoitns = [0.22665	0.45216	0.59526	-0.21777	-1.29612	-1.34484	-0.32403	-1.19616	-0.2445	-0.99195	-0.10563	-0.03651	0.03072	-0.17304	0.47913	0.56436	0.19569	-0.53937	-0.06246	-1.05588	-0.60489	-0.86127	-1.40427	0.59298	0.32502	-0.39462	0.07431	-0.30453	-0.05001	-0.33432	-0.78867	-1.36527	0.66441	0.12276	1.08648	0.2367	0.25968	0.49137	-0.32166	-0.69804	-0.94239	-1.49202	-1.27098	-0.60717	-0.23262	-1.05378	-0.23595	-1.48905	-1.5	0.6651	-0.92562	-0.26028	-1.02078	-1.03275	0.60843]

% then I continued the optimization using BBFGS which didn't work. So I
% went back to MOGA2 with more mutations. 
% zpoints = [ 0.40671	0.23865	0.52044	-0.15945	-1.25181	-1.42623	-0.10002	-0.80364	-0.72207	-1.10856	-0.18813	-0.03819	0.02271	-0.13383	0.45363	0.55449	0.02427	-0.69591	-0.14022	-1.10046	-0.56649	-1.00671	-1.26057	0.58848	0.2847	-0.4308	0.02832	-0.44166	-0.03693	-0.35625	-1.32675	-1.35687	0.62493	0.01125	1.09884	0.14595	0.26664	0.46449	-0.28599	-0.7326	-0.92916	-1.47702	-1.31127	-0.71127	-0.18636	-1.03698	-0.33162	-1.49571	-1.5	0.64182	-0.93726	-0.16728	-1.04703	-1.5	0.32868]

zpoints_v2(1,:) = zpoints(1:11);
zpoints_v2(2,:) = zpoints(12:22);
zpoints_v2(3,:) = zpoints(23:33);
zpoints_v2(4,:) = zpoints(34:44);
zpoints_v2(5,:) = zpoints(45:55);
doplot = 0;

[ydisplacment, cost,maxDisplacement] = FEALevelSet_2D_v7(xpoints,ypoints,zpoints_v2, doplot); % Call FEA for first time

obj = cost
cinequality = maxDisplacement-0.1;  % max Displacemen is positive. The max displacement must be less than 0.1 inch
cequality = ydisplacment+0.0623; % displacement will be negative, so add 0.0623 to make = 0
%feature('GetPid')





