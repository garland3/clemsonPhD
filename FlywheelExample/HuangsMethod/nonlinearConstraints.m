function [c,ceq] = nonlinearConstraints(x)
%c = ...     % Compute nonlinear inequalities at x.
% ceq = ...   % Compute nonlinear equalities at x.

[maxStress] = callflywheel2(x);
allowableStress = 45e6; % 45MPa

mass = callCaculateMass(x);
allowableMass = 75; % kg

c(1) = maxStress-allowableStress; % <= 0
c(2) = mass-allowableMass;  %<= 0

ceq = 0;



