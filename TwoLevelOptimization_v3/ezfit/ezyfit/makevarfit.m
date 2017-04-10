function makevarfit(f)
%MAKEVARFIT  Create variables from the parameters of a fit
%   MAKEVARFIT(F) creates (in the Matlab workspace) the variables that
%   contain the numerical values of the parameters from the fit F.
%
%   If you want the variables to be automatically created in the Matlab
%   workspace at each call of ezfit, showfit or selectfit, set the option
%   'automakevarfit = on' in fitparam.
%
%   If the input argument F is not specified, use the last fit.
%
%   Example:
%     Some sample data are fitted by a 2nd order polynom, and the
%     three variables 'a','b','c', which contain the numerical values
%     of the parameters, are created in the workspace:
%        plotsample('poly2');
%        f = showfit('a*x^2+b*x+c');
%        makevarfit(f);
%        whos
%
%   See also EZFIT, SHOWFIT, EDITCOEFF.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.00,  Date: 2006/02/15
%   This function is part of the EzyFit Toolbox

% History:
% 2006/02/15: v1.00, first version.


% if no input argument, use the last fit, which is stored
% in the variable lastfit in the 'base' workspace:
if nargin==0,
    if evalin('base','exist(''lastfit'',''var'')')
        f=evalin('base','lastfit');
    else
        errordlg('No existing fit coefficients. First fit a curve.',...
            'Edit Fit Coefficients','on');
        return;
    end;
end;

% create the variable f.param{i} containing the value f.m(i):
for i=1:length(f.param)
    assignin('base',f.param{i},f.m(i));
end;
