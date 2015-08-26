xpoints = 0:1:10; % X position of the control points

% From modeFrontier, meet y displacement constraint. Minimize the psi level
% set everywhere. Using simplex algorithm
% ypoints = [0.793504216	0.974776971	0.130891657	0.040603351	0.020277539	0.21006301	0.63956218	-1.370239927	-1.093348043	-1.499902779	-1.385100511]
% MOGA-II algorithm, 
% ypoints = [0.02235	0.41916	0.12546	0.17925	0.00396	-0.07017	-1.4994	-1.49982	-1.4817	-1.49994	-1.49997];
% NSGA2
% ypoints = [ 1.486871717	0.02486936658	0.5576749271	0.918654856	0.1815932342	0.3069125126	1.144076671	-0.4530789978	0.436152586	-0.1176224546	1.393037113]
% MOGA 2, with just 6 control points
% ypoints = [0.00345	1.17471	1.21566	-1.5	-1.49982	-1.49967];
% ---------------------------
% Version 3
% ypoints = [ 1.5	0.04083	0.11145	-0.25671	-0.62319	-0.41226]
%  version 3, where the first node was set to 0 so that the base would be
%  pla to help prevent warping.
% ypoints = [ 0	0.4086	0.22401	-0.29631	-0.81435	-0.25653]
% ---------------------------

% -------------------------
% Version 4
% ypoints = zeros(1,6) ; % displacement  is    -0.0418
% ypoints = ones(1,6) ; % displacement  is   -0.0828
% So the average displacement is  -0.0623 which is the new target
% displacement

% This is the one used in the journal paper. 
% Moga2 with DFM (pla at start)
%  ypoints = [0	0.75588	0.10422	-0.22041	-1.46037	-0.08451] % cost =  0.7155

% Moga 2 without DFM (first node does not have to be pla)
% ypoints = [ 0.80121	0.50874	0.62526	-0.76536	-0.76773	-0.35235] % cost =     0.7595


% ------------------
% Version 5
% Moga 2
 %ypoints = [0	1.0887	1.26678	-1.37712	-0.91218	-1.35783] % cost is 0.8368946958
 
 % With 10 control points
 % Moga2
 ypoints = [0	0.82323	0.29682	0.03189	0.00636	-0.00435	-0.33795	-0.72885	-0.1689	-1.05174	-0.39495]

doplot = 1;
[ydisplacment, cost]= FEALevelSet(xpoints,ypoints,doplot); % Call FEA for first time

