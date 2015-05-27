function GAOptimizationWrapper()
clear
clc
close all
rng default % for reproducibility

runningLocal = 0; % set equal to 1 if running on my local machine. Otherwise, set settings for Palmetto
runParallel = 1; % Set equal to 1 to run the algorithm in parallel
manualStartWorkers = 0; % set to 1 to manuall start works. 


nvar = 55;
LB = ones(1,55)*-3;   % Lower bound = -3
UB = ones(1,55)*3;   % Upper bound = 3

global xLast myf myc myceq count
xLast = []; % Last place computeall was called
myf = []; % Use for objective at xLast
myc = []; % Use for nonlinear inequality constraint
myceq = []; % Use for nonlinear equality constraint
count = 0; % count the number of times the analysis is called

fun = @objfun2; % the objective function, nested below
cfun = @constr2; % the constraint function, nested below


% Testing options
if(runningLocal ==1)
    gaoptions = gaoptimset('Generations',10,'Display','iter','PopulationSize',5, 'Elitecount', 2, 'TolCon',0.00623,'StallGenLimit',50);
else
    % actual values, I want to use. 
    gaoptions = gaoptimset('Generations',2000,'Display','diagnose','HybridFcn',@fmincon,  'PopulationSize',100, 'Elitecount', 5,'TolCon',0.00623, 'StallGenLimit',50);
end


% Run in Parallel
if(runParallel == 1)   
    
    if(manualStartWorkers ~=1)
        if(runningLocal ==1)
            numWorkersValue = 5;
        else
            numWorkersValue = 16;
        end

        c = parcluster;
        c.NumWorkers = numWorkersValue;  % where XXX is the number of works you want
        saveProfile(c);    
        matlabpool open
    end
    disp 'Number of workers'
    n = matlabpool('size') 
    gaoptions = gaoptimset(gaoptions,'UseParallel','always');

% Else, run in series
else
    % Special Settings for serial execution
    disp 'Serial execution of points'    
end

% General ga call by both
startTime = tic;
[x,fval, exitflag, output, population, scores] = ga(fun,nvar,[],[],[],[],LB,UB,cfun,gaoptions);
csvwrite('Xoutput.csv',x);
csvwrite('population.csv',population);
csvwrite('scores.csv',scores)
exitflag
x
fval
time_ga_sequential = toc(startTime);
fprintf(' optimization takes %g seconds.\n',time_ga_sequential);

if(runParallel == 1)
    if(manualStartWorkers ~=1)
     matlabpool close
    end
end

% diary off





end
