function y = objfun2_v2(x)
global xLast myf myc inequality


    if ~isequal(x,xLast) % Check if computation is necessary
        [myf,inequality] = FEALevelSetWrapperGA_v2(x);
        xLast = x;
    end
    % Now compute objective function
    y = myf;
end