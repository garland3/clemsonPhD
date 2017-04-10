function y=evalfit(f,x)
%EVALFIT  Evaluate a fit
%   Y = EVALFIT(F,X) evaluates the fit F for the values X.  F is a fit
%   structure, as obtained by EZFIT, SHOWFIT or SELECTFIT. Y is a vector
%   of the same length as X, with Y(i) = F(X(i)).
%
%   Example:
%     plotsample('power');
%     f = ezfit('power; log');
%     x = logspace(1,3,1000);
%     hold on; plot(x,evalfit(f,x),'r-'); hold off;
%
%   See also EZFIT, SHOWFIT, FITPARAM.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.01,  Date: 2006/09/04
%   This function is part of the EzyFit Toolbox

% History:
% 2006/02/16: v1.00, first version.
% 2006/09/04: v1.01, help improved.

m = f.m;
x__ref = x;
y = eval(eq2ml(f.eq, f.xvar));
if length(y)==1, % if the fit is a cste,
    y=y*ones(1,length(x));
end;