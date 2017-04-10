function dispeqfit(f,fp)
%DISPEQFIT   Display the equation of a fit.
%   DISPEQFIT(F) displays the equation of the fit F in the command window,
%   using the settings defined in FITPARAM. The fit structure F is
%   obtained from EZFIT. By default, DISPEQFIT is automatically called from
%   FIT when no output argument is specified.
%
%   Example:
%      plotsample('power')
%      f = ezfit('alpha/x^n');
%      dispeqfit(f);
%
%   See also FITPARAM, EZFIT, SHOWFIT, SHOWEQBOX

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.20,  Date: 2006/10/18
%   This function is part of the EzyFit Toolbox

% History:
% 2006/02/08: v1.00, first version.
% 2006/02/13: v1.10, compatible with 'y(x)=..' (free function name)
% 2006/10/18: v1.20, new argument fp

if nargin<2
    % loads the default fit parameters:
    try
        fp=fitparam;
    catch
        error('No fitparam file found.');
    end
end

streq = f.eq;
if strcmp(fp.eqreplacemode,'on'),
    for n=1:length(f.m),
        streq = strrep(streq, f.param{n}, num2str(f.m(n), fp.numberofdigit));
    end
    streq = strrep(streq,'+-','-');
    disp(['Equation: ' f.yvar '(' f.xvar ') = ' streq]);
else    
    disp(['Equation: ' f.yvar '(' f.xvar ') = ' streq]);
    for i=1:length(f.m),
        disp([ '     ' f.param{i} ' = ' num2str(f.m(i), fp.numberofdigit)]);
    end
end

lastline='';
switch lower(fp.corrcoefmode),
    case 'r', lastline=['R = ' num2str(f.r, fp.numberofdigit) '  '];
    case 'r2', lastline=['R^2 = ' num2str(f.r^2, fp.numberofdigit) '  '];
end
if strcmp(fp.linlogdisp,'on');
    lastline=[lastline '(' f.fitmode ')'];
end
if ~isempty(lastline)
    disp(['     ' lastline]);
end
