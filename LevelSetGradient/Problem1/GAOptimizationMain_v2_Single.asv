function GAOptimizationWrapper_v2()
% version 2 is the truss design similar to the classical Bendose paper
clear
clc
close all
rng default % for reproducibility

runningLocal = 0; % set equal to 1 if running on my local machine. Otherwise, set settings for Palmetto
runParallel = 1; % Set equal to 1 to run the algorithm in parallel
manualStartWorkers = 0; % set to 1 to manuall start works. 


nvar = 65;
LB = ones(1,nvar)*-1;   % Lower bound = -3
UB = ones(1,nvar)*1;   % Upper bound = 3

global xLast myf myc myceq count
xLast = []; % Last place computeall was called
myf = []; % Use for objective at xLast
myc = []; % Use for nonlinear inequality constraint
myceq = []; % Use for nonlinear equality constraint
count = 0; % count the number of times the analysis is called

fun = @objfun2_v2; % the objective function, nested below
cfun = @constr2_v2; % the constraint function, nested below
outPutFuct = @outPutFunction;

shrink = 0.7;
scale = 2;

% http://www.mathworks.com/help/gads/genetic-algorithm-options.html#f14223
% rows are each individual
% columns are the desgin vars

count2 = 1;
for temp  = 0:0.05:1
    initialP(count2,:) = ones(1,nvar)*temp;
    count2 = count2+1;
end

for temp2 = 1:nvar
    initialP(count2,temp)= (-1)^temp2;
end
count2 = count2+1;



% Testing options
if(runningLocal ==1)
    gaoptions = gaoptimset('Generations',10,'Display','iter','PopulationSize',5, 'Elitecount', 1,...
        'StallGenLimit',20,'OutputFcns',outPutFuct, ...
        'MutationFcn',{@mutationadaptfeasible, scale, shrink}, ...
        'InitialPopulation',initialP);
else
    % actual values, I want to use. 
    % population size = 5*23 where 23 is number of worker threads
    gaoptions = gaoptimset('Generations',2000,'Display','diagnose','HybridFcn',@fmincon,  ...
        'PopulationSize',100, 'Elitecount', 1,'StallGenLimit',20,'OutputFcns',outPutFuct,'MutationFcn',{@mutationadaptfeasible, scale, shrink});
end


% Run in Parallel
if(runParallel == 1)   
    
    if(manualStartWorkers ~=1)
        if(runningLocal ==1)
            numWorkersValue = 5;
        else
            numWorkersValue = 23; % minus 1 the number of cpus. Level one cpu for the main program. 
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
output.generations
output.message

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
