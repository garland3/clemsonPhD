xpoints = 0:2:10; % X position of the control points

% From modeFrontier, meet y displacement constraint. Minimize the psi level
% set everywhere. Using simplex algorithm
% ypoints = [0.793504216	0.974776971	0.130891657	0.040603351	0.020277539	0.21006301	0.63956218	-1.370239927	-1.093348043	-1.499902779	-1.385100511]

% MOGA-II algorithm, 
% ypoints = [0.02235	0.41916	0.12546	0.17925	0.00396	-0.07017	-1.4994	-1.49982	-1.4817	-1.49994	-1.49997];

% NSGA2
% ypoints = [ 1.486871717	0.02486936658	0.5576749271	0.918654856	0.1815932342	0.3069125126	1.144076671	-0.4530789978	0.436152586	-0.1176224546	1.393037113]


% MOGA 2, with just 6 control points
% ypoints = [0.00345	1.17471	1.21566	-1.5	-1.49982	-1.49967];
doplot = 0;

[ydisplacment, cost]= FEALevelSet(xpoints,ypoints,doplot); % Call FEA for first time

