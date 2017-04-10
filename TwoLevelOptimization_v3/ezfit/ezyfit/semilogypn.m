function [hp, hn] = semilogypn(x,y, specpos, specneg, varargin)
%SEMILOGYPN Semilogarithmic  plot for positive and negative data.
%   SEMILOGYPN(...) is the same as SEMILOGY(...), except that the Y data
%   may have positive and negative values. This is useful to detect
%   unexpected negative values in the Y data (causing repeated 'Negative
%   data ignored' warnings), or to visualize the magnitude of an
%   oscillating signal in log scale.
%
%   SEMILOGYPN(X,Y) plots X versus Y as logarithmic scales, where Y may
%   have positive and negative values. By default, negative values are
%   ploted with dashed lines and positive values with full lines. If X is
%   not specified, Y is plotted versus its index.
%
%   SEMILOGYPN(X,Y,LineSpecPos,LineSpecNeg) specifies LineSpecPos and
%   LineSpecNeg for the line types, marker symbols and colors for the
%   positive and negative values of Y. For instance,
%   SEMILOGYPN(X,Y,'ro-','b*') plots the positive Y with red line and 'o'
%   markers and the negative Y with blue stars '*'.
%   
%   SEMILOGYPN(X,Y,LineSpecPos,LineSpecNeg,'PropertyName',PropertyValue,...) 
%   sets property values for all lineseries graphics objects created by
%   SEMILOGYPN. See the line reference page for more information.
%
%   [HP, HN] = SEMILOGYPN(...) returns the handles to the two lineseries
%   graphics objects.
%
%   Example:
%      x = linspace(1,10,200);
%      y = sin(x*2)./x;
%      semilogypn(x,y,'r.','bo');
%      axis([1 10 1e-2 1]);
%
%      semilogypn(x,y,'b-','r:','LineWidth',2);
%
%      [hp, hn] = semilogypn(x,y,'k');
%      set(hp, 'LineWidth', 1);
%      set(hn, 'LineWidth', 2);
%
%   See also SEMILOGY, LOGLOG, LOGLOGPN.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.01,  Date: 2007/04/12
%   This function is part of the EzyFit Toolbox

% History:
% 2007/04/10: v1.00, first version.
% 2007/04/12: v1.01, do not draw lines between ranges of Y
%                    separated by ranges of opposite sign.

% error(nargchk(1,inf,nargin));

if nargin<2
    y = x;
    x = 1:length(y);
end

if nargin<3        % default LineSpecs:
    specpos = 'b-';
    specneg = 'b--';
elseif nargin<4    % LineSpec for Y<0 if the LineSpec for Y>0 is specified
    if findstr(specpos,'--')
        specneg=strrep(specpos,'--','-');
    else
        specneg=specpos;
        specneg=strrep(specneg,'-.','');
        specneg=strrep(specneg,':','');
        specneg=strrep(specneg,'-','');
        specneg=[specneg '--'];
    end
end

indpos = find(y>0);
indneg = find(y<0);

ypos = y;
ypos(indneg) = 0;

yneg = y;
yneg(indpos) = 0;

if isempty(indneg) && ~isempty(indpos)
    hp = semilogy(x, ypos, specpos, varargin{:});
    hn = [];
elseif ~isempty(indneg) && isempty(indpos)
    hp = [];
    hn = semilogy(x, -yneg, specneg, varargin{:});
else
    statehold = ishold;   % figure was 'held on' before?
    hp = semilogy(x, ypos, specpos, varargin{:});
    hold on
    hn = semilogy(x, -yneg, specneg, varargin{:});
    if statehold==false
        hold off
    end
end

if ~nargout
    clear hp hn
end
