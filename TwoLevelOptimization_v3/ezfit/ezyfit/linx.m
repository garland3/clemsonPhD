function linx
%LINX  Turn the X axis to LIN
%   LINX turns the X axis of the current plot to linear coordinates.
%
%   See also LINY, LOGX, LOGY.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.00,  Date: 2005/10/31
%   This function is part of the EzyFit Toolbox

set(gca,'XScale','lin');
