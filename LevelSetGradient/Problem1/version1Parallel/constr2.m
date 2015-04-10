function [c,ceq] = constr2(x)
global xLast myf myc myceq

        if ~isequal(x,xLast) % Check if computation is necessary
            [myf,myc,myceq] = FEALevelSetWrapperGA(x);
            xLast = x;
        end
        % Now compute constraint functions
        c = myc; % In this case, the computation is trivial
        ceq = myceq;
    end