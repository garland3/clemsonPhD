function docezyfit(filename)
%DOCEZYFIT  Documention for the Ezyfit toolbox
%   DOCEZYFIT displays the start page for the Ezyfit toolbox in the help
%   browser. In Matlab 7.3 and before (R2006), typing "doc ezyfit" had the
%   same result. However, since Matlab 7.4 (R2007a), the doc function
%   changed, and this feature is not available anymore, so you have to
%   use DOCEZYFIT instead.
%
%   DOCEZYFIT function_name  displays the documention of function_name
%   in the help browser. In Matlab 7.3, this is strictly equivalent to
%   DOC function_name.
%
%   Example:  DOCEZYFIT showfit

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.00,  Date: 2007/09/11
%   This function is part of the Ezyfit Toolbox

% History:
% 2007/09/11: v1.00, first version.

if nargin==0
    filename='ezyfit.html';
end


pathstr = fileparts(which('docezyfit'));
if strfind(filename,'.html')
    htmlfile = fullfile(pathstr, 'html', filename);
else
    htmlfile = fullfile(pathstr, 'html',  [filename '.html']);
end
web(htmlfile, '-helpbrowser');
