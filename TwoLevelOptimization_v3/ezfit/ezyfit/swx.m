function swx
%SWX  LIN<->LOG swap of the X axis
%   SWX swaps the X axis of the current plot between linear and log
%   coordinates.
%
%   See also SWY, SW.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.00,  Date: 2005/10/31
%   This function is part of the EzyFit Toolbox

if strcmp(get(gca,'XScale'),'log')
    linx;
else
    logx;
end
