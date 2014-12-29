function [c,ceq] = volumeconstraint(x)
%c = ...     % Compute nonlinear inequalities at x.
% ceq = ...   % Compute nonlinear equalities at x.

c = 0;
ceq = pi*x(1)^2*x(2)-10; % = 0