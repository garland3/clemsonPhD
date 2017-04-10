function checkupdate_ef(opt)
%CHECKUPDATE_EF  Check for update for Ezyfit
%   CHECKUPDATE_EF connects on the EzyFit webpage and checks for a new
%   version.
%
%   CHECKUPDATE_EF('dialog') does the same, but outputs the result in
%   a dialog box.
%
%   CHECKUPDATE_EF is automatically called at each Matlab restart if the
%   EzyFit menu is installed (see EFMENU).
%
%   See also ABOUT_EF, EFMENU. 

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.04,  Date: 2006/10/30
%   This function is part of the EzyFit Toolbox

% History:
% 2005/12/12: v1.00, first version.
% 2006/01/12: v1.01, opens a new web browser.
% 2006/04/10: v1.02, Toolbox renamed 'EzyFit'.
% 2006/10/03: v1.03, now checks the version number at the URL
%                    'www.fast.u-psud.fr/ezyfit/ezyfitversion.txt'
% 2006/10/30: v1.04, text in command window improved

verb=false;  % verbose mode

if nargin==0
    opt='command';
end

v1=ver('ezyfit');
curv=str2double(v1.Version);

[s,res]=urlread('http://www.fast.u-psud.fr/ezyfit/ezyfitversion.txt');
if ~res
    if strcmp(opt,'command')
        if verb, disp('Server unavailable.'); end;
        return
    else
        helpdlg('Server unavailable','Check for update');
        return
    end
end

newv=str2double(s);

if newv>curv
    if strcmp(opt,'command')
        disp(' ');
        disp(['   New: EzyFit ' s ' is now available online. <a href="matlab:web http://www.fast.u-psud.fr/ezyfit/html/ezyfit_releasenotes.html -browser">What''s new?</a>']);
        disp('   Click <a href="matlab:web http://www.fast.u-psud.fr/ezyfit -browser">here</a> to update.');
        disp(' ');
    else
        button = questdlg(['New: EzyFit ' s ' is now available online.'],'Check for update','Go online','What''s new?','Cancel','Go online');
        if strncmpi(button,'Go online',2),
            web http://www.fast.u-psud.fr/ezyfit -browser
        elseif strncmpi(button,'What',2),
            web http://www.fast.u-psud.fr/ezyfit/html/ezyfit_releasenotes.html -browser
        else
            return
        end
    end
elseif newv<curv
    if strcmp(opt,'command')
        if verb, disp('The online version is older than yours!'); end
    else
        helpdlg('The online version is older than yours!','Check for update');
    end
else
    if strcmp(opt,'command')
        if verb, disp('No new version available online'); end
    else
        helpdlg('No new version available online','Check for update');
    end
end
