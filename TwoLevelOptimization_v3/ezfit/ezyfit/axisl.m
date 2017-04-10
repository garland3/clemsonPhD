function axisl(opt)
%AXISL Enlarge the axis to the nearest power of 10.
%   AXISL enlarges the axis of the current plot to include the nearest
%   power of 10.
%
%   AXISL X or AXISL Y only includes the nearest power of 10 for the
%   horizontal or vertical axis.
%
%   See also AXISC, AXIS0, GRIDC.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.00,  Date: 2007/03/27
%   This function is part of the EzyFit Toolbox

% History:
% 2007/03/27: v1.00, first version.

if nargin==0
    opt='xy';
end

a=axis;

switch lower(opt)
    case 'x'
        axis([10^(floor(log10(a(1)))) 10^(ceil(log10(a(2)))) a(3) a(4)]);
    case 'y'
        axis([a(1) a(2) 10^(floor(log10(a(3)))) 10^(ceil(log10(a(4))))]);
    otherwise
        axis([10^(floor(log10(a(1)))) 10^(ceil(log10(a(2)))) 10^(floor(log10(a(3)))) 10^(ceil(log10(a(4))))]);
end

