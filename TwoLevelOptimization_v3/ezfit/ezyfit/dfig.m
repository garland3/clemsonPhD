function hh = dfig(h)
%DFIG  Create docked figure window
%   DFIG, by itself, creates a new docked figure window, and returns its
%   handle.
%
%   DFIG(H) makes H the current figure and docks it.  If Figure H does not
%   exist, and H is an integer, a new figure is created with handle H.
%
%   See also FIGURE.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.00,  Date: 2007/09/11


if nargin==0
    h=figure;    % create a new figure
else
    figure(h)    % makes H the current figure
end

set(h,'WindowStyle','docked');  % dock the figure

if nargout~=0   % returns the handle if requested
    hh=h;
end
