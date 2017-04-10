function swy
%SWY  LIN<->LOG swap of the Y axis
%   SWY swaps the Y axis of the current plot between linear and log
%   coordinates.
%
%   See also SWX, SW.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.00,  Date: 2005/10/31
%   This function is part of the EzyFit Toolbox

if strcmp(get(gca,'YScale'),'log')
    liny;
else
    logy;
end
