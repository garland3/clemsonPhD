function axis0(opt)
%AXIS0  Include the origin in the current axis
%   AXIS0 modifies the axis of the current plot to include the origin.
%   AXIS0 X or AXIS0 Y only includes the origin for the horizontal or
%   vertical axis.
%
%   See also AXISC, GRIDC.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.00,  Date: 2007/03/27
%   This function is part of the EzyFit Toolbox

% History:
% 2007/03/27: v1.00, first version.

if nargin==0,
    opt='xy';
end;

a=axis;

switch lower(opt)
    case'x',
        axis([min([0 a(1)]) max([0 a(2)]) a(3) a(4)]);
    case 'y'
        axis([a(1) a(2) min([0 a(3)]) max([0 a(4)])]);
    otherwise
        axis([min([0 a(1)]) max([0 a(2)]) min([0 a(3)]) max([0 a(4)])]);
end;

