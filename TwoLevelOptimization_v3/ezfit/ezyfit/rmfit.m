function rmfit(opt)
%RMFIT  Remove fits from the current figure
%   RMFIT removes all the fits (and their equation box) created by
%   SHOWFIT from the current figure.
%
%   RMFIT('first') removes only the first fit.
%   RMFIT('last') removes only the last fit (equivalent to UNDOFIT).
%   RMFIT(N) removes only the Nth fit.
%   RMFIT('boxonly') removes only the equation boxes.
%   RMFIT('getslope') removes all the equations displayed by GETSLOPE and
%   SHOWSLOPE.
%
%   Note that RMFIT does not remove the annotation objects (eg, lines).
%
%   See also SHOWFIT, UNDOFIT, GETSLOPE.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.10,  Date: 2006/02/17
%   This function is part of the EzyFit Toolbox

% History:
% 2005/10/31: v1.00, first version.
% 2005/11/03: v1.01, removes also the getslope textbox
% 2006/02/17: v1.10, updates the legend

% gr_dummy not defined for old versions
if verLessThan('matlab','8.4')
    eval('gr_dummy = 0;');
else
    gr_dummy = groot;
end

if exist('fitparam.m','file'),
    fp=fitparam;
else
    error('No fitparam file found.');
end

if nargin==0, opt=''; end

if isempty(get(gr_dummy,'CurrentFigure'))
    return
end

if isnumeric(opt)
    h=findall(gcf,'UserData','equationbox');
    h=h(end:-1:1);
    if length(h)>=opt, delete(h(opt)); end
    h=findobj(gcf,'UserData','fit');
    h=h(end:-1:1);
    if length(h)>=opt, delete(h(opt)); end
elseif ischar(opt)
    switch lower(opt)
        case 'first'
            h=findall(gcf,'UserData','equationbox');
            if ~isempty(h), delete(h(end)); end
            h=findobj(gcf,'UserData','fit');
            if ~isempty(h), delete(h(end)); end
        case 'last'
            h=findall(gcf,'UserData','equationbox');
            if ~isempty(h), delete(h(1)); end
            h=findobj(gcf,'UserData','fit');
            if ~isempty(h), delete(h(1)); end
        case 'boxonly'
            delete(findall(gcf,'UserData','equationbox'));
        case 'getslope'
            delete(findall(gcf,'UserData','getslopetextbox'));
            delete(findall(gcf,'UserData','showslopetextbox'));
        otherwise
            delete(findall(gcf,'UserData','equationbox'));
            delete(findobj(gcf,'UserData','fit'));
    end
end

% updates the legend:
if strcmp(fp.dispfitlegend,'on')
    legend off;  % refresh the legend
    legend show;
end
