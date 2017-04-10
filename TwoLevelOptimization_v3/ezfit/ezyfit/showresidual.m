function yres = showresidual(f,h)
%SHOWRESIDUAL   Show the residuals of a fit
%   SHOWRESIDUAL(F) displays the residuals for the fit F in a new figure,
%   ie: plots the difference between the data and the fit.
%
%   SHOWRESIDUAL is also available from the item 'Show Fit Residuals' in
%   the EzyFit menu (see EFMENU).
%
%   If F is not specified, plots the residual for the last fit.
%
%   YRES = SHOWRESIDUAL(...) returns the residuals (YRES = YDATA - YFIT).
%
%   SHOWRESIDUAL(F, H) plots the residual in the figure H.
%
%   See also SHOWFIT, EZFIT, EFMENU.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.00,  Date: 2006/02/16
%   This function is part of the EzyFit Toolbox

% History:
% 2006/02/16: v1.00, first version.


% if no input argument, use the last fit, which is stored
% in the variable lastfit in the 'base' workspace:
if nargin==0,
    if evalin('base','exist(''lastfit'',''var'')')
        f=evalin('base','lastfit');
    else
        errordlg('No existing fit. First fit a curve.',...
            'Show fit residual','on');
        return;
    end;
end;

yfit = evalfit(f, f.x);
yres = f.y - yfit;

if exist('h','var'),
    figure(h);
    hold on;
else
    figure;
end;

hres = plot(f.x, yres, 'o');
title(['Fit Residuals for ' strrep(f.name,'^','\^')]);
ylabel([f.yvar ' - ' f.yvar '_{fit}']);

if isfield(f,'hdata')
    if ishandle(f.hdata),
        set(hres, 'Color', get(f.hdata, 'Color'));
        set(hres, 'LineStyle', get(f.hdata, 'LineStyle'));
        set(hres, 'LineWidth', get(f.hdata, 'LineWidth'));
        set(hres, 'Marker', get(f.hdata, 'Marker'));
        set(hres, 'MarkerSize', get(f.hdata, 'MarkerSize'));
        p=get(f.hdata); pp=get(p.Parent); ppl=get(pp.XLabel);
        xlabel(ppl.String); % The Xlabel of the residual plot is given
        % by the XLabel of the original plot.
        set(gca, 'XScale', pp.XScale);
    end;
end;

gridc x;
axisc y;
hold off;

if nargout==0,
    clear yres;
end;
