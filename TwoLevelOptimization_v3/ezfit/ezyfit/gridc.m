function gridc(opt)
%GRIDC  Centered cross grid.
%   GRIDC shows central cross axes.
%   GRIDC X or GRIDC Y only shows the horizontal or vertical axe.
%
%   See also AXISC, AXIS0.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.01,  Date: 2006/02/16
%   This function is part of the EzyFit Toolbox

% History:
% 2005/11/10: v1.00, first version.
% 2006/02/16: v1.01, option 'x' or 'y'.

if nargin==0,
    opt='xy';
end;

a=axis;
hold on;
if findstr(lower(opt),'x');
    plot([a(1) a(2)],[0 0],'k:');
end;
if findstr(lower(opt),'y');
    plot([0 0],[a(3) a(4)],'k:');
end;
hold off;
