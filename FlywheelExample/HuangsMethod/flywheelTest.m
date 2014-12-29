function [mass, energy, volume, maxStress] = flywheelTest()
close all;
clear;clc;


innerR =0.02 % m,  inner shoudl be 0.02
outerR = 0.2 % m,  out should be 0.2

% Optimization 5.12 in Huang, results
% hpoints = [innerR,0.132205485913699;
%             0.04,0.0532472762500307;
%              0.08,0.021923194;
%               0.12,0.048987798;
%              outerR,0.052049772]         
% 
% VolFractpoints = [innerR,0.540397065;
%             0.04, 0.97383401  ;
%              0.08,  0.245426554;
%               0.12,   0.394520906;
%              outerR,  0.217602125]
         
% hpoints = [innerR,1;
%             0.04,1;
%              0.08,1;
%               0.12,1;
%              outerR,1]         
% 
% VolFractpoints = [innerR,1;
%             0.04, 1  ;
%              0.08, 1;
%               0.12, 1;
%              outerR, 1]



% HeightControlPoints[0]	HeightControlPoints[1]	HeightControlPoints[2]	HeightControlPoints[3]	HeightControlPoints[4]	VolFractSnControlPoints[0]	VolFractSnControlPoints[1]	VolFractSnControlPoints[2]	VolFractSnControlPoints[3]	VolFractSnControlPoints[4]	Energy	Mass	MaxStress	obj_stress	energy_constraint	mass_constraint	<ERROR>	<VIRTUAL>	<FEASIBLE>
% 0.2	0.02	0.02	0.02	0.2	0.0	0.0	0.0	0.0	0.0	150422.7669386106	31.81835198845218	3.2558850118075658E7	3.2558850118075658E7	150422.7669386106	31.81835198845218	false	false	true
hpoints = [innerR,0.2;
            0.04,0.02;
             0.08,0.02;
              0.12,0.02;
             outerR,0.2]         

VolFractpoints = [innerR,0.0;
            0.04, 0.0  ;
             0.08, 0.0;
              0.12, 0;
             outerR, 0]
%          
[mass, energy, volume, maxStress] = flywheel(hpoints, VolFractpoints)