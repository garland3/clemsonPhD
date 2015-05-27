function [c,ceq] = constr2list7(x)
global xLast myf  inequality

        if ~isequal(x,xLast) % Check if computation is necessary
            [myf,inequality] = FEALevelSetWrapperGA_v2(x);
            xLast = x;
        end
        % Now compute constraint functions
        c = inequality-0.7; % In this case, the computation is trivial
        ceq = []; % for the truss design problem, there are no equality constraints
    end