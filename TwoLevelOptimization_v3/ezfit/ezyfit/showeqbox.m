function hboxeq=showeqbox(f,varargin)
%SHOWEQBOX   Show the equation box of a fit.
%   H = SHOWEQBOX(F) shows the equation box of the fit F, using the settings
%   defined in FITPARAM. The fit structure F is obtained from EZFIT. H is a
%   handle to the equation box.
%
%   Note that SHOWEQBOX interprets the greek symbols in the fitting equation
%   (latex syntax). The '\' latex symbol must be omitted.
%
%   SHOWEQBOX is automatically called from SHOWFIT or SELECTFIT when the
%   option 'dispeqboxmode' is set to 'on' in the fitparam.m file.
%
%   Example:
%      plotsample('power')
%      f = ezfit('alpha/x^n');
%      showeqbox(f);
%
%   See also FITPARAM, EZFIT, SHOWFIT, DISPEQFIT.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.14,  Date: 2010/07/09
%   This function is part of the EzyFit Toolbox

% History:
% 2006/02/08: v1.00, first version.
% 2006/02/14: v1.10, compatible with 'y(x)=..' (free function name).
%                    Capital greek letters accepted. Greek letters as
%                    subscript accepted.
% 2006/09/06: v1.11, bug fixed when the handle f.hdata is invalid
% 2006/10/18: v1.12, bug fixed for fitcolors. accepts the argument fp
% 2007/07/24: v1.13, bug fixed when no 2nd parameter given
% 2010/07/09: v1.14, now the font size and font name are the same as the
%                    axes of the figure.

% gr_dummy not defined for old versions
if verLessThan('matlab','8.4')
    eval('gr_dummy = 0;');
else
    gr_dummy = groot;
end

if nargin>1   % new v1.13
    if isstruct(varargin{1})
        fp=varargin{1};
    end
end

if ~exist('fp','var')
    % loads the default fit parameters:
    try
        fp=fitparam;
    catch
        error('No fitparam file found.');
    end
end

if isfield(f,'eq'), % for normal fits (no interpolation):
    streq = [f.yvar '(' f.xvar ') = ' f.eq];
    
    if strcmp(fp.eqreplacemode,'on'),
        for n=1:length(f.m),
            streq = strrep(streq, f.param{n}, num2str(f.m(n), fp.numberofdigit));
        end
        if length(streq)>fp.maxlengtheq,
            streq=[streq(1:fp.maxlengtheq) '...'];
        end
        streq=greekize(streq);
        streq={streq};
    else
        % truncates the equation string if too long:
        if length(streq)>fp.maxlengtheq,
            streq=[streq(1:fp.maxlengtheq) '...'];
        end
        streq={greekize(streq)};
        for n=1:length(f.m),
            strm = [greekize(f.param{n}) ' = ' num2str(f.m(n), fp.numberofdigit)];
            streq = {streq{:} strm};
        end
    end
    
    lastline='';
    switch lower(fp.corrcoefmode)
        case 'r', lastline=['R = ' num2str(f.r, fp.numberofdigit) '  '];
        case 'r2', lastline=['R^2 = ' num2str(f.r^2, fp.numberofdigit) '  '];
    end
    if strcmp(fp.linlogdisp,'on')
        lastline=[lastline '(' f.fitmode ')'];
    end
    if ~isempty(lastline)
        streq = {streq{:} lastline};
    end
    
else
    streq=f.name; % for interpolations
end



% number of fit already present in the figure:
numann=length(findall(gcf,'UserData','equationbox'));
% position of the new annotation textbox, to avoid overlapping:
position = fp.boxlocation+numann*[0.01 -0.01 0 0];

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

if any(strncmpi(varargin,'transparent',5)) % option used with 'selectfit', in order not to hide the data during the selection
    bgcolor='none';
else
    bgcolor='white';
end

% changed v1.14
hboxeq=annotation(...
    'textbox',position,...
    'BackgroundColor',bgcolor,...
    'Color',fitcolor,...
    'EdgeColor',fitcolor,...
    'FitBoxToText','on',...
    'UserData','equationbox',...
    'String',streq,...
    'FontSize',get(gr_dummy,'DefaultAxesFontSize'),...
    'FontName',get(gr_dummy,'DefaultAxesFontName'));

