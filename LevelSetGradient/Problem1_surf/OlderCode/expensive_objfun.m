function f = expensive_objfun(x)
%EXPENSIVE_OBJFUN An expensive objective function used in optimparfor example.

%   Copyright 2007-2012 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2012/08/21 01:09:29 $

% Simulate an expensive function by performing an expensive computation
eig(magic(300));
% Evaluate objective function
f = exp(x(1)) * (4*x(3)^2 + 2*x(4)^2 + 4*x(1)*x(2) + 2*x(2) + 1);
