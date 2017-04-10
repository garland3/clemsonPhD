function logy
%LOGY  Turn the Y axis to LOG
%   LOGY turns the Y axis of the current plot to log coordinates.
%
%   See also LOGX, LINX, LINY.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.00,  Date: 2005/10/31
%   This function is part of the EzyFit Toolbox

set(gca,'YScale','log');
