function liny
%LINY  Turn the Y axis to LIN
%   LINY turns the Y axis of the current plot to linear coordinates.
%
%   See also LINX, LOGX, LOGY.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.00,  Date: 2005/10/31
%   This function is part of the EzyFit Toolbox

set(gca,'YScale','lin');
