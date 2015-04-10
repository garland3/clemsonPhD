function gaParallel()

xLast = []; % Last place computeall was called
myf = []; % Use for objective at xLast
myc = []; % Use for nonlinear inequality constraint
myceq = []; % Use for nonlinear equality constraint

fun = @objfun; % the objective function, nested below
cfun = @constr; % the constraint function, nested below



%matlabpool open
% n = matlabpool('size')
%if max(size(gcp)) == 0 % parallel pool needed
%    matlabpool % create the parallel pool
%end
nvar = 4;
gaoptions = gaoptimset('Generations',15,'Display','iter');


rng default % to get the same evaluations as the previous run
gaAvailable = true;
if gaAvailable
    gaoptions = gaoptimset(gaoptions,'UseParallel','always');
    startTime = tic;
    [x,f,eflag,outpt] = ga(fun,nvar,[],[],[],[],[],[],cfun ,gaoptions);
    x
    f
    eflag
    outpt
    time_ga_parallel = toc(startTime);
    fprintf('Parallel GA optimization takes %g seconds.\n',time_ga_parallel);
end

% matlabpool close

function y = objfun(x)
    if ~isequal(x,xLast) % Check if computation is necessary
        [myf,myc,myceq] = computeall(x);
        xLast = x;
    end
    % Now compute objective function
    y = myf + 20*(x(3) - x(4)^2)^2 + 5*(1 - x(4))^2;
end

function [c,ceq] = constr(x)
    if ~isequal(x,xLast) % Check if computation is necessary
        [myf,myc,myceq] = computeall(x);
        xLast = x;
    end
    % Now compute constraint functions
    c = myc; % In this case, the computation is trivial
    ceq = myceq;
end


end