function [n,a]=getslope(varargin)
%GETSLOPE  Slope of the current line.
%   GETSLOPE displays the equation of the current line of the figure. Use
%   the menu 'Insert > Line' to draw a line first. GETSLOPE allows for
%   rough curve fitting "by eye". You may also use the shortcut Ctrl+G if
%   the Ezyfit menu has been installed (see EFMENU).
%
%   Depending of the axis types, the equation of the line will be:
%        Y = N*X+A         for X linear and Y linear
%        Y = A*X^N         for X log and Y log
%        Y = A*EXP(N*X)    for X linear and Y log
%        Y = A+N*LOG(X)    for X log and Y linear (LOG = natural logarithm)
%
%   GETSLOPE('Property1',...) specifies the display mode:
%        'figure'   output the result in the figure (by default)
%        'command'  output the result in the command window
%        'nodisplay'  no output
%        'slope'    displays only the slope N (by default).
%        'equation' displays the full equation (e.g., N*X+A)
%
%   [N, A] = GETSLOPE(...) also returns the parameters N and A of the
%   equation.
%
%   See also SHOWSLOPE, GETLINEINFO, PLOTSAMPLE, RMFIT, EFMENU.


%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.41,  Date: 2011/10/04
%   This function is part of the EzyFit Toolbox


% History:
% 2005/02/17: First version.
% 2005/02/23: uses the current active line
% 2005/02/24: bug 'nodisp' fixed
% 2005/05/21: standard error messages
% 2005/05/24: option 'draw' suppressed; works in linear or log coordinates.
% 2005/05/26: Bug fixed (zero slope in semilogx mode)
% 2005/06/16: Displays the fit parameters with the correct number of digits.
% 2005/06/22: v1.23, bug fixed in semilogx mode; 3 digits.
% 2005/07/27: v1.24, correctly displays the negative values (without '+').
% 2005/11/02: v1.30, works also when one extremity of the line is selected.
%                    the equation is now in an annotation box, that follows
%                    the line.
% 2006/02/08: v1.31, use the number of digits defined in fitparam
% 2006/03/08: v1.32, option 'dialog' added, when called from efmenu
% 2006/05/09: v1.33, bug fixed for exponential fits (annotation box too
%                    small)
% 2007/09/14: v1.40, the procedure for the detection of the line has been
%                    rewritten for ML 7.4-7.5; new options 'command' and
%                    'figure'
% 2011/10/04: v1.41, bug fixed for matlab version >=7.10 (thanks Alex !!)

% gr_dummy not defined for old versions
if verLessThan('matlab','8.4')
    eval('gr_dummy = 0;');
else
    gr_dummy = groot;
end

% loads the default fit parameters:
try
    fp=fitparam;
catch
    fp.numberofdigit=3;
end


if isempty(get(gr_dummy,'CurrentFigure'))
    error('No figure.');
end

if isempty(legend)  % new v1.40
    islegend = false;
else
    islegend = true;
    legend off
end

strv=version;
v=str2double(strv(1:3));

if v>=7.4 || v == 7.1                 % changed v1.41 (thanks Alex!!)
    % get infos about the current line
    if ~ishandle(gco)
        if any(strncmpi(varargin,'dialog',4))
            errordlg('First select a line (Menu Insert > Line)','GetSlope');
            return
        end
        error('First select a line (Menu Insert > Line).');
    end
    
    h=gco; % currently active line
    o=get(h);
    if ~isfield(o,'X')
        if any(strncmpi(varargin,'dialog',4))
            errordlg('First select a line (Menu Insert > Line)','GetSlope');
            return
        end
        error('First select a line (Menu Insert > Line).');
    end
    % this is the coord of the line, normalized to the window (between 0 and 1)
    xw1=o.X(1); xw2=o.X(2);
    yw1=o.Y(1); yw2=o.Y(2);
    
else   % ------------------- FOR MATLAB <=7.3
    % get infos about the current line
    if ishandle(gco),
        % gco is the currently active line (but it does not work
        % if the line has just being drawn!)
        if ~isfield(get(gco),'XData')
            % gco is a handle, but there is no XData field in it
            if (gco==gcf),
                % nothing is selected...
                if any(strncmpi(varargin,'dialog',4))
                    errordlg('First select a line (Menu Insert > Line)','GetSlope');
                    return;
                end
                error('First select a line (Menu Insert > Line).');
            else
                % this is the case when the line was just being drawn.
                % (h(1) is the XAxis, h(2) and h(3) are the first and second points,
                % and h(4) is the line itself).
                h=findall(gcf,'Type','line');
                if length(h)<4,
                    if any(strncmpi(varargin,'dialog',4))
                        errordlg('First select a line (Menu Insert > Line)','GetSlope');
                        return;
                    end
                    error('First select a line (Menu Insert > Line).');
                else
                    h=h(4);
                end
            end
        else
            % the object has a XData field in it
            h=gco;
            o=get(gco);
            if ~(length(o.XData)==2),
                if (length(o.XData)==1),
                    % One point of the line was being moved.
                    % go back to its parent and take the 3rd child:
                    % here you will find the line itself.
                    o=get(o.Parent);
                    h=o.Children(3);
                else
                    if any(strncmpi(varargin,'dialog',4))
                        errordlg('First select a line (Menu Insert > Line)','GetSlope');
                        return;
                    end
                    error('First select a line (Menu Insert > Line).');
                end
            end
        end
    else
        if any(strncmpi(varargin,'dialog',5))
            errordlg('First select a line (Menu Insert > Line)','GetSlope');
            return;
        end
        error('First select a line (Menu Insert > Line).');
    end
    o=get(h);
    % this is the coord of the line, normalized to the window (between 0 and 1)
    xw1=o.XData(1); xw2=o.XData(2);
    yw1=o.YData(1); yw2=o.YData(2);
    
end


icg=get(gca); % get infos about the current axes

% computes the coordinates of the line
% normalized to the axes of the figure (between 0 and 1):
%   icg.Position(1:2) = coordinates (x,y) relative to the window
%   icg.Position(3:4) = width and height of the windows
xr1=(xw1-icg.Position(1))/(icg.Position(3));
xr2=(xw2-icg.Position(1))/(icg.Position(3));


yr1=(yw1-icg.Position(2))/(icg.Position(4));
yr2=(yw2-icg.Position(2))/(icg.Position(4));

xb1=icg.XLim(1);
xb2=icg.XLim(2);
if strcmpi(icg.YDir, 'reverse')
    yb2=icg.YLim(1);  yb1=icg.YLim(2);
else
    yb1=icg.YLim(1);  yb2=icg.YLim(2);
end

if strncmpi(icg.XScale,'log',3) && strncmpi(icg.YScale,'log',3)
    % Fit by a power law, y=a*x^n.
    
    % computes the coordinates of the line in 'physical' units (those of
    % the axes):
    x1=xb1*(xb2/xb1)^xr1;
    x2=xb1*(xb2/xb1)^xr2;
    
    y1=yb1*(yb2/yb1)^yr1;
    y2=yb1*(yb2/yb1)^yr2;
    
    n=log(y2/y1)/log(x2/x1);
    a=y1/(x1^n);
    
    if any(strncmpi(varargin,'equation',2))
        str=[num2str(a, fp.numberofdigit) ' x^{' num2str(n, fp.numberofdigit) '}'];
    else
        str=num2str(n, fp.numberofdigit);
    end
    
elseif strncmpi(icg.XScale,'linear',3) && strncmpi(icg.YScale,'linear',3)
    % Fit by an affine law, y=nx+a.
    
    x1=xb1+(xb2-xb1)*xr1;
    x2=xb1+(xb2-xb1)*xr2;
    
    y1=yb1+(yb2-yb1)*yr1;
    y2=yb1+(yb2-yb1)*yr2;
    
    n=(y2-y1)/(x2-x1);
    a=y1-n*x1;
    
    if any(strncmpi(varargin,'equation',2))
        if a>0 % changed 27/07/2005, v1.24
            str=[num2str(n, fp.numberofdigit) ' x + ' num2str(a, fp.numberofdigit)];
        else
            str=[num2str(n, fp.numberofdigit) ' x ' num2str(a, fp.numberofdigit)];
        end
    else
        str=num2str(n, fp.numberofdigit);
    end
    
elseif strncmpi(icg.XScale,'linear',3) && strncmpi(icg.YScale,'log',3)
    % Fit by an exponential law, y=a*exp(nx).
    
    x1=xb1+(xb2-xb1)*xr1;
    x2=xb1+(xb2-xb1)*xr2;
    
    y1=yb1*(yb2/yb1)^yr1;
    y2=yb1*(yb2/yb1)^yr2;
    
    n=log(y2/y1)/(x2-x1);
    a=y1/exp(n*x1);
    
    if any(strncmpi(varargin,'equation',2))
        str=[num2str(a, fp.numberofdigit) ' e^{' num2str(n, fp.numberofdigit) ' x}'];
    else
        str=num2str(n, fp.numberofdigit);
    end
    
elseif strncmpi(icg.XScale,'log',3) && strncmpi(icg.YScale,'linear',3)
    % Fit by a logarithmic law, y=a+n*log(x).
    
    x1=xb1*(xb2/xb1)^xr1;
    x2=xb1*(xb2/xb1)^xr2;
    
    y1=yb1+(yb2-yb1)*yr1;
    y2=yb1+(yb2-yb1)*yr2;
    
    n=(y2-y1)/log(x2/x1);
    if n~=0
        a=y1-n*log(x1);
    else
        a=NaN;
    end
    if any(strncmpi(varargin,'equation',2))
        str=[num2str(a, fp.numberofdigit) ' ln (' num2str(n, fp.numberofdigit) ' x)'];
    else
        str=num2str(n, fp.numberofdigit);
    end
end


% if the current line has a textbox previously attached to it, kill it:
if ishandle(get(h,'UserData'))
    delete(get(h,'UserData'));
end

if ~any(strncmpi(varargin,'nodisplay',3))
    if any(strncmpi(varargin,'command',3))
        disp(str)
    end
    
    if (~any(strncmpi(varargin,'figure',3)) && ~any(strncmpi(varargin,'command',3))) ...
            || any(strncmpi(varargin,'figure',3))
        textlocation=[(xw1+xw2)/2 (yw1+yw2)/2-0.075-sign(n)*0.025 0.4 0.1];  % changed v1.33
        htext=annotation('textbox',textlocation,'LineStyle','none','UserData','getslopetextbox','String',str);
        % says to the line that the textbox htext is attached to it:
        set(h,'UserData',htext);
    end
end

if islegend
    legend show   % put back the legend (new v1.40)
end

if nargout==0
    clear n
end
