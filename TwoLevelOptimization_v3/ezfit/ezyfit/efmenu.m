function efmenu(opt)
%EFMENU   Ezyfit menu
%   EFMENU adds or refreshes the Ezyfit menu for the current figure and
%   for all new figures.
%   EFMENU OFF removes the menu from all the figures and for all new
%   figures.
%
%   Careful: when the Ezyfit menu is active, all saved figure files (.FIG)
%   include the Ezyfit menu. This generates a problem when those figure
%   files are opened in a Matlab system which does not have the Ezyfit
%   toolbox installed. In order to remove the Ezyfit menu from the figure
%   file, use the function REMOVE_EFMENU_FIG.
%
%   If you want to always have the EzyFit menu in your figures, type
%   EFMENU INSTALL. This will create or update the 'startup.m' file in the
%   user directory of your Matlab installation. In addition, at
%   each Matlab restart, this will check the last version of the EzyFit
%   toolbox on the web (see CHECKUPDATE_EF).
%
%   In order to cancel the effect of EFMENU INSTALL, you need to edit the
%   'startup.m' file and to remove manually the 2 lines of code related to
%   the Ezyfit toolbox.
%
%   See also PLOTSAMPLE, SHOWFIT, EZFIT, UNDOFIT, RMFIT, EDITFIT,
%   LOADFIT, CHECKUPDATE_EF, REMOVE_EFMENU_FIG.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.90,  Date: 2016/04/28
%   This function is part of the EzyFit Toolbox

% History:
% 2005/12/06: v1.00, first version.
% 2005/12/16: v1.01, works with all new figures.
% 2006/01/12: v1.10, User defined fits improved.
% 2006/01/19: v1.20, User defined fits improved again.
% 2006/01/30: v1.21, bug fixed with selectfit. New fit fgauss.
% 2006/02/08: v1.30, use a new mat-file for user defined fits.
%                    polynomial orders as a submenu. New 'Edit Fit Coeff'.
%                    Check for update when launched from startup.m
% 2006/02/16: v1.31, item 'Show Fit Residual' added
% 2006/02/27: v1.40, option 'install' added.
% 2006/03/07: v1.41, option 'install' improved (append to startup.m), and
%                    bug fixed for check of the ml version
% 2006/04/10: v1.42, toobox renamed 'EzyFit'
% 2006/04/24: v1.43, bugs in option 'install' fixed (did not correctly
%                    recognize if the toolbox was already installed)
% 2006/09/06: v1.44, submenu for 'PlotSample'
% 2007/09/12: v1.45, new menu 'Refresh EzyFit menu'
% 2007/09/14: v1.46, bug fixed with 'install off' (case sensitive); options
%                    modified for getslope and showslope
% 2007/09/18: v1.50, new submenus for Gaussian fits
% 2008/03/06: v1.51, Ctrl-Z removed (conflicts)
% 2008/03/14: v1.60, selectfit removed for Ezyfit >= 2.30 (replaced by
%                    showfit in "Brushed" mode)
% 2009/02/04: v1.61, help text 'uninstall' added
% 2010/07/07: v1.70, help text 'remove_efmenu_fig' added
% 2012/06/28: v1.71, check update removed
% 2012/07/02: v1.80, change working directory, for Windows 7 compatibility
% 2012/08/01: v1.81, change for under MacOS X, Matlab 2012 (thanks J. Bagrow)
% 2014/03/25: v1.82, bug fixed for Windows (C:\)
% 2014/10/24: v1.90, bug fixed for Matlab 2014b for new groot behavior
%                   (thanks to Olaf Bousché)
% 2016/04/28: v1.91, changed text in startup file



% check the matlab version:
% if str2double(version('-release'))<14,
%     error('EzyFit requires Matlab 7 (R14) or higher.');
% end


% gr_dummy not defined for old versions (new v1.90)
if verLessThan('matlab','8.4')
    eval('gr_dummy = 0;');
else
    gr_dummy = groot;
end

% if efmenu is called from the file startup.m, then calls
% the Check Update for the Ezyfit Toolbox:
st=dbstack;
if length(st)>=2
    if strcmp(st(2).name,'startup'),
        checkupdate_ef;
    end
end


if nargin==0
    opt='on';
end

% new v1.30:
userfit=loadfit('user');

switch lower(opt)
    case 'off'
        set(gr_dummy,'DefaultFigureCreateFcn','');
        delete(findobj('Label','EzyFit')); % delete the menu from all the figures.
    case 'install'
        install_ef;
    case 'on'
        set(gr_dummy,'DefaultFigureCreateFcn','efmenu');
        if ~isempty(get(gr_dummy,'CurrentFigure'))
            delete(findobj(gcf,'Label','EzyFit')); % delete the menu only from the current figure
            h=uimenu('Label','EzyFit');
            uimenu(h,'Label','Undo fit','Callback','undofit'); % Ctrl-Z removed march 7, 2008
            uimenu(h,'Label','Remove all fits','Callback','rmfit');
            
            %submenu showfit and selectfit:
            hm(1)=uimenu(h,'Separator','on','Label','Show Fit');
            menuname{1}='showfit';
            %hm(2)=uimenu(h,'Label','Select Fit');
            %menuname{2}='selectfit';
            for i=1:1, % only "showfit" is visible for EzyFit >=2.30
                uimenu(hm(i),'Label','Constant','Callback',['lastfit=' menuname{i} '(''cste'');']);
                hlin=uimenu(hm(i),'Label','Linear');
                uimenu(hlin,'Label','a*x','Callback',['lastfit=' menuname{i} '(''linear'');']);
                uimenu(hlin,'Label','a*x+b','Callback',['lastfit=' menuname{i} '(''affine'');']);
                uimenu(hlin,'Label','a*(x-b)','Callback',['lastfit=' menuname{i} '(''affineshift'');']);
                hpol=uimenu(hm(i),'Label','Polynomial');
                for or=1:6,
                    uimenu(hpol,'Label',['Order ' num2str(or)],...
                        'Callback',['lastfit=' menuname{i} '(''poly' num2str(or) ''');']);
                end
                uimenu(hpol,'Label','Other...','Separator','on','Callback',['lastfit=' menuname{i} '(''poly'');']);
                hexp=uimenu(hm(i),'Label','Exponential');
                uimenu(hexp,'Label','a*exp(b*x)','Callback',['lastfit=' menuname{i} '(''exp'');']);
                uimenu(hexp,'Label','a*exp(b*x)+c','Callback',['lastfit=' menuname{i} '(''expc'');']);
                uimenu(hexp,'Label','a*exp(x/b)','Callback',['lastfit=' menuname{i} '(''expdiv'');']);
                uimenu(hexp,'Label','a*(1-exp(-x/b))','Callback',['lastfit=' menuname{i} '(''explim'');']);
                hlog=uimenu(hm(i),'Label','Logarithmic');
                uimenu(hlog,'Label','a*log(b*x)','Callback',['lastfit=' menuname{i} '(''log'');']);
                uimenu(hlog,'Label','a*log(x)+b','Callback',['lastfit=' menuname{i} '(''logc'');']);
                hpow=uimenu(hm(i),'Label','Power');
                uimenu(hpow,'Label','a*x^n','Callback',['lastfit=' menuname{i} '(''power'');']);
                uimenu(hpow,'Label','a*x^n+c','Callback',['lastfit=' menuname{i} '(''powerc'');']);
                uimenu(hpow,'Label','a*(x+b)^n','Callback',['lastfit=' menuname{i} '(''powershift'');']);
                hgauss=uimenu(hm(i),'Label','Gaussian');
                uimenu(hgauss,'Label','Centered Normalized Gaussian','Callback',['lastfit=' menuname{i} '(''cngauss'');']);
                uimenu(hgauss,'Label','Centered Non-normalized Gaussian','Callback',['lastfit=' menuname{i} '(''cfgauss'');']);
                uimenu(hgauss,'Label','Non-centered Normalized Gaussian','Callback',['lastfit=' menuname{i} '(''ngauss'');']);
                uimenu(hgauss,'Label','Non-centered Non-normalized Gaussian','Callback',['lastfit=' menuname{i} '(''fgauss'');']);
                hosc=uimenu(hm(i),'Label','Oscillations');
                uimenu(hosc,'Label','a*sin(b*x)','Callback',['lastfit=' menuname{i} '(''sin'');']);
                uimenu(hosc,'Label','a*cos(b*x)','Callback',['lastfit=' menuname{i} '(''cos'');']);
                uimenu(hosc,'Label','a*sin(b*x)+c','Callback',['lastfit=' menuname{i} '(''sinc'');']);
                uimenu(hosc,'Label','a*cos(b*x)+c','Callback',['lastfit=' menuname{i} '(''cosc'');']);
                uimenu(hosc,'Label','a*sin(b*x+phi)','Callback',['lastfit=' menuname{i} '(''sinphi'');']);
                uimenu(hosc,'Label','a*sin(b*x+phi)+c','Callback',['lastfit=' menuname{i} '(''sinphic'');']);
                uimenu(hm(i),'Label','Nearest','Separator','on','Callback',[menuname{i} '(''nearest'');']);
                uimenu(hm(i),'Label','Spline','Callback',[menuname{i} '(''spline'');']);
                uimenu(hm(i),'Label','Cubic/pchip','Callback',[menuname{i} '(''pchip'');']);
                
                for j=1:length(userfit)
                    if j==1,
                        uimenu(hm(i),'Label',['#' num2str(j) ': ' userfit(j).name],...
                            'Separator','on','Callback',['lastfit=' menuname{i} '(''' userfit(j).eq ''');']);
                    else
                        uimenu(hm(i),'Label',['#' num2str(j) ': ' userfit(j).name],...
                            'Callback',['lastfit=' menuname{i} '(''' userfit(j).eq ''');']);
                    end
                end
                uimenu(hm(i),'Label','Other...','Separator','on','Callback',[menuname{i} ';']);
            end
            
            %submenu edit fit:
            hef=uimenu(h,'Label','Edit User Fit');
            for j=1:length(userfit)
                uimenu(hef,'Label',['#' num2str(j) ': ' userfit(j).name],'Callback',['editfit(' num2str(j) ')']);
            end
            uimenu(hef,'Label','New User Fit...','Separator','on','Callback','editfit');
            uimenu(hef,'Label','Reset','Separator','on','Callback',...
                'if strcmp(questdlg(''Are you sure you want to reset the user defined fits?'',''Reset User Fits'',''No''),''Yes''), delete(''userfit.mat''); efmenu; end');
            
            
            % remainder of the menu:
            uimenu(h,'Label','Refresh EzyFit menu','Callback','efmenu');
            uimenu(h,'Label','Edit Fit Coefficients','Callback','editcoeff');
            uimenu(h,'Label','Show Fit Residuals','Callback','showresidual');
            
            uimenu(h,'Label','Get Slope','Separator','on','Callback','getslope(''equation'',''dialog'')','Accelerator','G');
            uimenu(h,'Label','Show Slope...','Callback','showslope([],''dialog'')');
            
            hps=uimenu(h,'Label','Plot Sample','Separator','on');
            uimenu(hps,'Label','Linear','Callback','plotsample linear');
            uimenu(hps,'Label','Constant','Callback','plotsample cste');
            uimenu(hps,'Label','3 polynomial curves','Callback','plotsample poly2');
            uimenu(hps,'Label','Oscillations','Callback','plotsample osc');
            uimenu(hps,'Label','Damped oscillations','Callback','plotsample damposc');
            uimenu(hps,'Label','Exponential decay','Callback','plotsample exp');
            uimenu(hps,'Label','Histogram with one peak','Callback','plotsample hist');
            uimenu(hps,'Label','Histogram with two peaks','Callback','plotsample hist2');
            uimenu(hps,'Label','Power law','Callback','plotsample power');
            uimenu(hps,'Label','Power law with a cutoff','Callback','plotsample powco');
            
            hqc=uimenu(h,'Label','Axis','Separator','on');
            uimenu(hqc,'Label','Swap X-axis lin <-> log','Callback','swx');
            uimenu(hqc,'Label','Swap Y-axis lin <-> log','Callback','swy');
            uimenu(hqc,'Label','Swap both X- and Y- axis','Callback','sw','Accelerator','B');
            uimenu(hqc,'Label','Include the origin','Callback','axis0');
            uimenu(hqc,'Label','Center the origin','Callback','axisc');
            uimenu(hqc,'Label','Enlarge to nearest power of 10','Callback','axisl');
            
            uimenu(h,'Label','Edit Default Settings...','Separator','on','Callback','edit fitparam');
            uimenu(h,'Label','Help EzyFit','Callback','docezyfit');
            uimenu(h,'Label','Check for update','Callback','checkupdate_ef(''dialog'')');
            uimenu(h,'Label','EzyFit Home Page','Callback','web http://www.fast.u-psud.fr/ezyfit -browser');
            uimenu(h,'Label','About EzyFit','Callback','about_ef(''dialog'')');
        end
    otherwise
        error('Unknown option');
end


% -----------------------------------------
% Subfunction 'efmenu install'


function install_ef

% This location no longer works for Windows 7:
%sufile = fullfile(matlabroot,'toolbox','local','startup.m');

sufolder = strrep(userpath,';','');
sufolder = strrep(sufolder,':',''); % new v1.81, thanks J. Bagrow
sufolder = strrep(sufolder,'C\','C:\'); % new 1.82, was buggy for windows
sufile = fullfile(sufolder,'startup.m');

ft = [];
if exist(sufile,'file')
    try
        ft = textread(sufile,'%s');
    catch
        ft = [];
    end
end
if ~isempty(ft)
    for num = 1:length(ft)
        if strfind(ft{num},'efmenu'),
            disp('The EzyFit menu is already installed in your startup.m file:');
            disp(sufile);
            return
        end
    end
    %copyfile(sufile,fullfile(matlabroot,'toolbox','local','startup_previous.m'));
    copyfile(sufile,fullfile(sufolder,'startup_previous.m'));
    fid=fopen(sufile,'a');
    fprintf(fid,'%s\n',' ');
    fprintf(fid,'%s\n',['%   These lines have been added by ''efmenu install'' (' datestr(now) '):']);
    fprintf(fid,'%s\n','efmenu;   % Includes the EzyFit menu for all new figure.');
    fprintf(fid,'%s\n\n','fprintf('' <a href="matlab:docezyfit">EzyFit</a> - A free curve fitting toolbox for Matlab.\n\n'');');
    fclose(fid);
    
    disp('The EzyFit menu has been correctly installed in your startup.m file.');
    disp('To get started, select <a href="matlab:docezyfit">EzyFit</a>.');
    
else % creates a new startup.m file and re-calls install_ef:
    fid=fopen(sufile,'w');
    fprintf(fid,'%s\n','%STARTUP   Startup file');
    fprintf(fid,'%s\n','%   This file is executed when MATLAB starts up.');
    fprintf(fid,'%s\n',' ');
    fclose(fid);
    install_ef;
end
