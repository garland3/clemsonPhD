% Copyright Anthony Garland 2015
function  [history,searchdir] = optimize(w1,w2,iterationNum)
LASTN = maxNumCompThreads(8)

% convert the command line arguments from strings to doubles
w1 = str2double(w1)
w2 = str2double(w2)



% flywheel optimizatoin using fmincon
% do equation 5.15 in Huang's disertation

% Set up shared variables with OUTFUN
history.x = [];
history.fval = [];
searchdir = [];


% Inputs are the bezier curve control points. 
% Half are used for the profile and half are used for the volume fraction
% composition. 
numberInputs = 16;

half = numberInputs/2;
A = zeros(numberInputs,numberInputs);
b = zeros(numberInputs,1);
Aeq = [];
beq = [];
lb(1:half) = 0.02; % lower bounds for the profile control points
ub(1:half) = 0.1; % upper bounds for the profile control points
lb(half+1:numberInputs) = 0; % lower bounds for the vol frac control points
ub(half+1:numberInputs) = 1; % upper bounds for the vol frac  control points
x0(1:half) = (0.1-0.02)/2; % start in the middle of profile control points bounds
x0(half+1:numberInputs) = 0.5; % start in the middle of vol frac control points bounds
options =   optimoptions('fmincon','MaxFunEvals',1000,'TolCon',1,'TolFun',100,'TolX',0.01,'UseParallel','Always','Display','iter'); %,'OutputFcn',@outfun,'Display','iter','Algorithm','active-set');
results = zeros(11,numberInputs);

message = sprintf('Optimize where w1 = %4.2f and w2 = %4.2f',w1,w2);
disp(message);

% run fmincon
[results,fval] = fmincon(@(x) Objective(x,w1,w2),x0,A,b,Aeq,beq,lb,ub,@nonlinearConstraints, options);
filename= strcat('results',num2str(iterationNum), '.csv'); % generate the filename
csvwrite(filename,results) % write the results to a file. 
results


end

