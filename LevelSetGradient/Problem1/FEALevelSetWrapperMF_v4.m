% xpoints = 0:1:10; % X position of the control points
[xpoints,ypoints] = meshgrid(0:2:10,0:0.5:1);

zpoints_v2 = ones(3,6)



% These are the design varriables. The optimizer controls these  z heights
%zpoints = 0:17
%zpoints = peaks(xpoints,ypoints);

% Version 6 of the FEA, not incorporates a 2D level set function with
% control points. 
% These results don't make since
% cost =  0.7075, displacement =  -0.0598
% zpoints = [0	0.21477	-0.48579	0.3645	0.14652	-1.49853	0	1.5	1.49118	-1.5	-1.49997	-1.37463	0	-0.58404	-1.5	-1.5	-1.5	-1.49526]
%zpoints = [0	-0.0582	-1.49454	-1.25691	0.71418	-0.17763	0	0.23877	1.48683	-1.49901	-1.38378	-1.48146	0	0.24381	-1.5	0.96102	-1.5	-1.49943]

%  cost = 0.7241
% zpoints = [0	-0.22998	-1.5	-1.26456	0.69342	-0.95175	0	0.91005	1.48431	-0.80916	-1.42056	-1.47864	0	0.15123	-1.5	0.85854	-1.4427	-1.49997]
zpoints = [ 0	-0.01761	-1.5	-1.5	-0.12588	-1.0698	0	0.76488	1.48953	-1.36128	-1.45662	-1.5	0	-0.34083	-1.49922	0.89559	-1.4997	-1.2336]

zpoints_v2(1,:) = zpoints(1:6)
zpoints_v2(2,:) = zpoints(7:12)
zpoints_v2(3,:) = zpoints(13:18)
doplot = 1;
[ydisplacment, cost,maxDisplacement]= FEALevelSet_2D(xpoints,ypoints,zpoints_v2, doplot) % Call FEA for first time

