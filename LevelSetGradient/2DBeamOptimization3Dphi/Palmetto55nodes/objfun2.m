function y = objfun2(x)
global xLast myf myc myceq


    if ~isequal(x,xLast) % Check if computation is necessary
        [myf,myc,myceq] = FEALevelSetWrapperGA(x);
        xLast = x;
    end
    % Now compute objective function
    y = myf;
end