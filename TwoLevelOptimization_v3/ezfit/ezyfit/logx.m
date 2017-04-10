function logx
%LOGX  Turn the X axis to LOG
%   LOGX turns the X axis of the current plot to log coordinates.
%
%   See also LINX, LINY, LOGY.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.00,  Date: 2005/10/31
%   This function is part of the EzyFit Toolbox

set(gca,'XScale','log');
