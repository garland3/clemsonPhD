function [sigmaVM,yDisp] = beamTheory(x,y)
%This function gives the beam theory solution for von mises stresses at a 
%point with the given xy coordinates. It includes shear force.
%I have checked this against a rework of the problem in mathCAD and it
%gives the same results.
L = 10; % in length
t = 0.1; % in thickness
h = 1; % in height
I = 1/12*t*h^3; %in^4 second moment of area
rho = 0.28; %lb/in^3 density
E = 29000000;
FappliedLoad = 0.25;
Mapplied = FappliedLoad*(L-x); %moment at the point from the applied load
Mbody = rho*t*h*(L-x)^2/2; %moment from body force
Mtot = Mapplied+Mbody; %total moment from body and applied
yNuet = y-h/2; %the y position relative to the nuetral axis.
sigmaX = Mtot*yNuet/I; %standard M*y/I from statics.
tauApplied = FappliedLoad/(t*h);%V=F/A Shear stress for the applied load.
tauBody = (L-x)*rho; %(L-x)*rho*t*h/(t*h) Shear stress from weight.
tau = tauApplied+tauBody;
sigmaVM = 1/2*sqrt(2*sigmaX^2+6*tau^2);
dispBody = t*h*rho*x^2/(24*E*I)*(4*L*x-x^2-6*L^2);%y Displacement due to body force.
dispApplied = FappliedLoad*x^2/(6*E*I)*(x-3*L);%y Displacement due to Applied force.
yDisp = dispBody+dispApplied;%total y Displacement
end

