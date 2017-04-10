function fout = showfit(varargin)
%SHOWFIT   Show the fit for the active curve.
%   SHOWFIT displays the fit specified in a dialog box for the active
%   curve. Use UNDOFIT or RMFIT to remove the fit. The display settings are
%   defined in the file FITPARAM.
%
%   By default, the first curve in the active figure is used (see FITPARAM
%   to change this default behavior). To fit another curve, select it
%   before calling SHOWFIT. If some data are selected by the "Data
%   Brushing" tool (only for Matlab >= 7.6), only those data are fitted.
%
%   SHOWFIT(FUN) specifies the fitting string FUN. FUN may be either a
%   default fit (eg, 'exp', 'power'), a user-defined fit (see EDITFIT),
%   or directly a fit equation (eg, 'c+a*exp(-x/x0)'). See EZFIT for the
%   syntax of FUN. FUN may also be any valid interpolation method string
%   (eg 'spline', 'cubic'... excepted 'linear'). See INTERP1 for the valid
%   interpolation methods.
%
%   SHOWFIT(F) displays the fit F defined from the fit structure returned
%   by EZFIT. See EZFIT for the definition of F.
%
%   SHOWFIT(..., 'PropertyName','PropertyValue',...) specifies the 
%   properties of the fit (e.g., fit colors, line width etc.). See the
%   default values in the file FITPARAM.
%
%   F = SHOWFIT(...) also returns the fit structure F. F has the same
%   content as the fit structure returned by EZFIT (see EZFIT for details),
%   and also contains a handle to the equation box and to the curve.
%
%   Note that if the option 'lin' or 'log' is not specified in the string
%   FUN, SHOWFIT checks the Y-scale to know if Y or LOG(Y) has to be
%   fitted.
%
%   Examples:
%      type 'plotsample' and follow the instructions.
%
%      plotsample power
%      showfit('c*x^n; log');
%
%      plotsample hist
%      f = ezfit('gauss');
%      showfit(f,'fitcolor','red','fitlinewidth',2);
%
%      plotsample poly2
%      f = showfit('z(v) = poly3','dispfitlegend','on');
%      editcoeff(f);
%
%   See also EZFIT, UNDOFIT, RMFIT, PLOTSAMPLE, PICKDATA, INTERP1.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 2.11,  Date: 2012/07/08
%   This function is part of the EzyFit Toolbox

% History:
% 2005/05/21: v1.00, first version.
% 2005/07/27: v1.01, cosmetics.
% 2005/10/07: v1.02, also displays the parameters in the plot.
% 2005/10/17: v1.10, uses the first curve of the active figure if no curve
%                    is active.
% 2005/10/18: v1.11, no display if no output arg.
% 2005/10/31: v1.20, also displays R. text is displayed in an annotation
%                    box. now uses pickdata. also accepts interpolation.
% 2005/11/03: v1.30, use the file 'fitparam.m' for the fit default values.
% 2005/11/05: v1.31, also displays the lin/log mode
% 2005/12/06: v1.40, opens a dialog box if the fitting function string is
%                    not specified.
% 2006/01/28: v1.41, first check if a curve is present.
% 2006/01/31: v1.42, bug fixed for R^2
% 2006/02/08: v2.00, new syntax. now works with the fit structure F.
% 2006/02/10: v2.01, displays the legend.
% 2006/02/16: v2.02, use evalfit.
% 2006/03/08: v2.03, idem v1.40; bug fixed when using showfit(f) where f
%                    is obtained from fit(x,y,...) without plot.
% 2006/09/04: v2.04, '$' and '!' accepted instead of ';' in FUN
% 2006/09/06: v2.05, bug fixed when F is a structure of a FIT originating
%                    from a deleted curve (f.hdata is an invalid handle)
% 2006/10/18: v2.10, accepts additional input parameters to chnage the
%                    default settings.
% 2012/07/08: v2.11, hold on/off now correctly (thanks Alec de Zegher)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit parameters:  (new v2.10)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loads the default fit parameters:
try
    fp=fitparam;
catch
    error('No fitparam file found.');
end


% change the default values of the fit parameters according to the
% additional input arguments:
for nopt=1:(nargin-1)
    if any(strcmp(varargin{nopt},fieldnames(fp))) % if the option string is one of the fit parameter
        fp.(varargin{nopt}) = varargin{nopt+1};
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input arguments:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==0,
    % opens a dialog box with the last fit equation (saved in prevfit.mat):
    efroot=fileparts(mfilename('fullpath'));   % directory where the ezyfit toolbox is installed
    prevfitfile=[efroot filesep 'prevfit.mat'];
    if exist(prevfitfile,'file'),
        load(prevfitfile);
    else
        fun='m1 + m2 * x; m1 = 1; m2 = 1';
    end
    fun=inputdlg('Enter a default fit, a user-defined fit or directly a fit equation',...
        'General Curve Fit Definition',1,{fun});
    if ~isempty(fun)
        fun=fun{1};
        fun=strrep(fun,'$',';');  % new v2.04
        fun=strrep(fun,'!',';');
        save(prevfitfile,'fun');
        % if the option 'lin' or 'log' is not specified, use the Y-scale of
        % the current figure:
        if (isempty(findstr(strrep(fun,' ',''),';lin'))) && (isempty(findstr(strrep(fun,' ',''),';log'))),
            mode=get(gca,'YScale');
            fun=strrep([fun ';' mode(1:3)],';;',';');
        end
        f=ezfit(fun);
    else
        return;
    end
elseif nargin>=1,
    if isstruct(varargin{1}), % if argument is a structure
        f=varargin{1};
    elseif ischar(varargin{1}), % if argument is a string
        fun=varargin{1};
        fun=strrep(fun,'$',';');  % new v2.04
        fun=strrep(fun,'!',';');
        if sum(strcmp(fun,{'nearest','spline','pchip','cubic','v5cubic'})),
            f.name=fun;
            [f.x, f.y, f.hdata] = pickdata(fp);
        else  % if the string is not an interpolation
            if strcmp(fun,'poly'), % new v2.03
                str_ord=inputdlg('Order of the polynomial fit','Polynomial order',1,{'2'});
                if ~isempty(str_ord),
                    fun=['poly' str_ord{1}];
                else
                    return
                end
            end
            % if the option 'lin' or 'log' is not specified, use the Y-scale of
            % the current figure:
            if (isempty(findstr(strrep(fun,' ',''),';lin'))) && (isempty(findstr(strrep(fun,' ',''),';log'))),
                mode=get(gca,'YScale');
                fun=strrep([fun ';' mode(1:3)],';;',';');
            end
            f=ezfit(fun,varargin{2:end});
        end
    else
        error('Use showfit with a string or a fit structure.');
    end
end



% output fit structure: idem as f, but with additional information.
fout=f;

% determines the X points used to evaluate the fit
switch fp.extrapol
    case {'none','data'},
        xminfit=min(f.x);
        xmaxfit=max(f.x);
    case 'fig',
        icg=get(gca);
        xminfit=icg.XLim(1);
        xmaxfit=icg.XLim(2);
end
if strcmp(get(gca,'XScale'),'log'),
    xfit=logspace(log10(xminfit),log10(xmaxfit),fp.npt);
else
    xfit=linspace(xminfit,xmaxfit,fp.npt);
end


% evaluate the fit / interpolation:
switch f.name
    case {'nearest','spline','pchip','cubic','v5cubic'}
        yfit = interp1(f.x, f.y, xfit, f.name);
    otherwise   % real fit (ie with a fitting function)  
        yfit = evalfit(f, xfit);
end

        
% determines the fit color:
if ischar(fp.fitcolor) || length(fp.fitcolor)==3
    fitcolor=fp.fitcolor; % fixed color
else
    fitcolor=[0 0 0]; % default color if no data in the figure
    if isfield(f,'hdata')
        if ishandle(f.hdata)
            co=get(f.hdata); % object properties of the data
            if isfield(co,'Color'),
                fitcolor=max(0,min(1,co.Color*fp.fitcolor)); % color indexed from that of the data
            end
        end
    end
end

% plots the fitted curve (modified v2.11)
holdState = ishold;
hold on;
fout.hfit = plot(xfit,yfit,'-','Color',fitcolor,'LineStyle',fp.fitlinestyle,'LineWidth',fp.fitlinewidth);
% Return to original hold state
if ~holdState
    hold off
end

% this tag is useful to delete the fit later
set(fout.hfit,'UserData', 'fit'); 

% give a 'DisplayName' to the fit:
if length(fout.name)>10,
    fitname='fit';
else
    fitname=greekize(fout.name);
end
setautodisplayname; % give a name to all the data in the figure (if they don't have one yet)
if isfield(f,'hdata'),
    if ishandle(f.hdata)   % new v2.05, changed v2.10
        try
            dataname = get(f.hdata,'DisplayName');
        catch
            dataname = 'data';
        end
        fitname=[fitname ' (' dataname ')'];
    end
end
set(fout.hfit,'DisplayName',fitname);


% makes the equation box:
if strcmp(fp.dispeqboxmode,'on')
    fout.heqbox = showeqbox(f,fp);
end

% updates the legend:
if strcmp(fp.dispfitlegend,'on')
    legend off;  % refresh the legend
    legend show;
end

% opens the Array Editor with the fit coefficients:
if (strcmp(fp.editcoeffmode,'on') && (isfield(fout,'m'))),
    editcoeff(fout);
end

% displays the fit equation in the command window if no output argument:
if ~nargout
    if isfield(fout,'m'),
        if strcmp(fp.dispeqmode,'on') % new v2.30
            dispeqfit(fout,fp);
        end
    end
    clear fout;
end
