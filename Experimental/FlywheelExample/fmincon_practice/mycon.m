function [c,ceq] = mycon(x)
%c = ...     % Compute nonlinear inequalities at x.
% ceq = ...   % Compute nonlinear equalities at x.

c = x(1)^2+x(2)^2+50; % <= 0
ceq = 0;

end

