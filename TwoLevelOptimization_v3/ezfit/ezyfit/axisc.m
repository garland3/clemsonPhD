function axisc(opt)
%AXISC  Centered axis
%   AXISC centers the axis of the current plot.
%   AXISC X or AXISC Y only centers the horizontal or vertical axe.
%
%   See also AXIS0, GRIDC.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.02,  Date: 2006/02/28
%   This function is part of the EzyFit Toolbox

% History:
% 2005/11/10: v1.00, first version.
% 2006/02/16: v1.01, option 'x' or 'y'.
% 2006/02/28: v1.02, bug fixed for option 'xy'.

if nargin==0,
    opt='xy';
end;

a=axis;
xbound = max([abs(a(1)) abs(a(2))]);
ybound = max([abs(a(3)) abs(a(4))]);

switch lower(opt)
    case'x',
        axis([-xbound xbound a(3) a(4)]);
    case 'y'
        axis([a(1) a(2) -ybound ybound]);
    otherwise
        axis([-xbound xbound -ybound ybound]);
end;

