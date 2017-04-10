function hl=showslope(n,varargin)
%SHOWSLOPE  Draw a line with fixed slope
%   SHOWSLOPE does not work under recent versions of Matlab - sorry!
%
%   SHOWSLOPE(N) drags a line in the current figure with a fixed slope N
%   (does not work with docked windows; undock your window first!)
%
%   Depending of the axis types, the 'slope' N means:
%        Y = N*X+A         for X linear and Y linear
%        Y = A*X^N         for X log and Y log
%        Y = A*EXP(N*X)    for X linear and Y log
%        Y = A+N*LOG(X)    for X log and Y linear (LOG = natural logarithm)
%
%   The value N is displayed close to the line. N may also be the string of
%   any valid Matlab expression (eg, 'pi', '22/7', ...), which will be
%   displayed close to the line. SHOWSLOPE(N,'nolabel') does not display
%   the value N in the figure.
%
%   H = SHOWSLOPE(..) also returns a handle to the line.
%
%   By default, the line is an 'annotation object', ie: it is attached to
%   the window and not to the figure axes. As a consequence, the line may
%   be further moved, or may be used for GETSLOPE. However, resizing
%   the window may in some case shift the line, and turning the axes from
%   lin to log turns a power law to a log law... To avoid these problems,
%   specify SHOWSLOPE(N,'fix'), which makes the line a real plot object
%   attached to the figure axes.
%
%   See also GETSLOPE, GETLINEINFO, PLOTSAMPLE.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.50,  Date: 2014/10/24
%   This function is part of the EzyFit Toolbox

% History:
% 2005/02/17, v1.00, first version.
% 2005/02/22, v1.01, cosmetic changes.
% 2005/04/20: v1.02, id.
% 2005/07/27: v1.03, id.
% 2005/09/03: v1.04, check arg changed.
% 2005/10/31: v1.05, uses myginput.
% 2005/11/02: v1.06, suppresses the extreme points; N can be a string.
% 2005/11/03: v1.10, annotation line used. options 'd' and 's' suppressed.
%                    works in lin/log coords.
% 2005/12/06: v1.20, opens a dialog box if the slope is not specified.
% 2005/12/12: v1.21, interactive properties of the figure suspended during
%                    the drawing.
% 2006/01/27: v1.30, uses 'SelectionType' instead of 'WindowButtonDownFcn',
%                    which caused frequent bugs.
% 2006/02/16: v1.40, saves the slope in the file 'prevslope.mat'
% 2007/09/14: v1.41, error msg when window is docked; new option 'dialog'
% 2014/10/24: v1.50, bug fixed for Matlab 2014b for new groot behavior
%                   (thanks to Olaf Bousché)

% error(nargchk(0,2,nargin));

% gr_dummy not defined for old versions (v1.50)
if verLessThan('matlab','8.4')
    eval('gr_dummy = 0;');
else
    gr_dummy = groot;
end

if isempty(get(gr_dummy,'CurrentFigure'))
    error('No figure.');
end

% new v1.41
if strcmpi(get(gcf,'WindowStyle'),'Docked')
    if any(strncmpi(varargin,'dialog',4))
        errordlg('Sorry, ShowSlope does not work on docked windows. Please undock your window first.','ShowSlope');
        return
    end
    error('Sorry, showslope does not work on docked windows. Please undock your window first.');
end

if nargin==0 || isempty(n)
    efroot=fileparts(mfilename('fullpath'));   % directory where the ezyfit toolbox is installed
    prevslopefile=[efroot filesep 'prevslope.mat'];
    if exist(prevslopefile,'file')
        load(prevslopefile);
    else
        n=1;
    end
    n=inputdlg('Enter the slope','Showslope',1,{num2str(n)});
    if ~isempty(n)
        n=n{1};
        save(prevslopefile,'n');
    else
        return;
    end
end


if ischar(n), strn=n; n=eval(n);
else strn=num2str(n); end

icg=get(gca); % get infos about the current axes

% initial settings:
fig=gcf;
initial_pointer=get(fig,'Pointer');
uistate=uisuspend(fig);

[xv,yv,but]=myginput(1,'crosshair');

set(fig,'Pointer','crosshair');

if (but==1),
    spt0 = get(gr_dummy, 'PointerLocation'); % pointer coordinates in pixel, from the screen corner
    wls=get(fig,'Position'); % position of the window in the screen, in pixels
    xw1=(spt0(1)-wls(1))/wls(3); % 1st point, normalized to the window (between 0 and 1)
    yw1=(spt0(2)-wls(2))/wls(4);
    
    hold on;
    
    set(fig,'SelectionType','extend');
    while strcmp(get(fig,'SelectionType'),'extend'),
        if get(gr_dummy,'PointerWindow')==fig, % if the pointer is in the window
            sptcur = get(gr_dummy, 'PointerLocation');    % pointer coordinates in pixel, from the screen corner
            wls=get(fig,'Position');               % position of the window in the screen, in pixels
            xw2=(sptcur(1)-wls(1))/wls(3);         % abscisse of the current point, normalized to the window
            xr2=(xw2-icg.Position(1))/(icg.Position(3));      % abscisse of the current point, normalized to axes
            
            if isequal(icg.XScale,'log') && isequal(icg.YScale,'log'),
                % Show a power law, y=a*x^n.
                
                xv(2)=icg.XLim(1)*(icg.XLim(2)/icg.XLim(1))^xr2;         % in 'physical' units
                yv(2)=yv(1)*(xv(2)/xv(1))^n;                             % in 'physical' units
                yr2=log(yv(2)/icg.YLim(1))/log(icg.YLim(2)/icg.YLim(1)); % normalized to the axes
                
            elseif isequal(icg.XScale,'linear') && isequal(icg.YScale,'linear'),
                % Show an affine law, y=nx+a.
                
                xv(2)=icg.XLim(1)+(icg.XLim(2)-icg.XLim(1))*xr2;
                yv(2)=yv(1)+n*(xv(2)-xv(1));
                yr2=(yv(2)-icg.YLim(1))/(icg.YLim(2)-icg.YLim(1));
                
            elseif isequal(icg.XScale,'linear') && isequal(icg.YScale,'log'),
                % Show an exponential law, y=a*exp(nx).
                
                xv(2)=icg.XLim(1)+(icg.XLim(2)-icg.XLim(1))*xr2;
                yv(2)=yv(1)*exp(n*(xv(2)-xv(1)));
                yr2=log(yv(2)/icg.YLim(1))/log(icg.YLim(2)/icg.YLim(1));
                
            elseif isequal(icg.XScale,'log') && isequal(icg.YScale,'linear'),
                % Show a logarithmic law, y=a+n*log(x).
                
                xv(2)=icg.XLim(1)*(icg.XLim(2)/icg.XLim(1))^xr2; % in 'physical' units
                yv(2)=yv(1)+n*log(xv(2)/xv(1));
                yr2=(yv(2)-icg.YLim(1))/(icg.YLim(2)-icg.YLim(1));
                
            end
            
            % normalization to the window (between 0 and 1):
            if n~=0
                yw2=icg.Position(2)+icg.Position(4)*yr2;
            else
                yw2=yw1;
            end
            hl=findall(fig,'UserData','showslopedashedline');
            if ~isempty(hl),
                set(hl,'X',[xw1 xw2]);
                set(hl,'Y',[yw1 yw2]);
            else
                annotation('line',[xw1 xw2],[yw1 yw2],...
                    'LineStyle',':','UserData','showslopedashedline');
            end
            %drawnow;
        end
    end
    
    if any(strcmpi(varargin,'fix'))
        delete(hl);
        x=linspace(xv(1),xv(2),200);
        if     isequal(icg.XScale,'log') && isequal(icg.YScale,'log'),       y=yv(1)*(x/xv(1)).^n;
        elseif isequal(icg.XScale,'linear') && isequal(icg.YScale,'linear'), y=yv(1)+n*(x-xv(1));
        elseif isequal(icg.XScale,'linear') && isequal(icg.YScale,'log'),    y=yv(1)*exp(n*(x-xv(1)));
        elseif isequal(icg.XScale,'log') && isequal(icg.YScale,'linear'),    y=yv(1)+n*log(x/xv(1));
        end
        hl=plot(x,y,'k-','UserData','showslopeline');
    else
        set(hl,'UserData','showslopeline');
        set(hl,'LineStyle','-');
    end
    %drawnow;
    
    if ~any(strncmpi(varargin,'nolabel',5))
        textlocation=[(xw1+xw2)/2 (yw1+yw2)/2-0.075-sign(n)*0.025 0.1 0.1];
        annotation('textbox',textlocation,'LineStyle','none','UserData','showslopetextbox','String',strn);
    end
end
hold off;

% restores all:
set(fig,'Pointer',initial_pointer);
uirestore(uistate);

if nargout==0
    clear hl
end
