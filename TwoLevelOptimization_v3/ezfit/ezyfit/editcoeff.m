function c=editcoeff(f)
%EDITCOEFF  Edit the coefficients of a fit
%   EDITCOEFF opens the Matlab's Array Editor to edit the coefficients
%   of the last fit (which is stored in the variable 'lastfit'). This may
%   be useful to copy/paste the fit coefficients into a spreadsheet
%   (see the setting 'coeffarray' in FITPARAM to change the row/column
%   display).
%
%   EDITCOEFF is also available from the item 'Edit Fit Coefficients' in
%   the EzyFit menu, and can be automatically called after each fit (see
%   the setting 'editcoeffmode' in FITPARAM).
%
%   EDITCOEFF(F) opens the Matlab's Array Editor to edit the coefficients
%   of the fit F. (F is a fit structure as returned by EZFIT or SHOWFIT). 
%
%   C = EDITCOEFF(...) returns a 2xN cell array, containing the
%   coefficients names and values.
%
%   Example:
%     plotsample('power');
%     f = ezfit('power');
%     editcoeff(f);
%
%   See also EZFIT, SHOWFIT, EFMENU, FITPARAM, MAKEVARFIT, OPENVAR.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.00,  Date: 2006/02/09
%   This function is part of the EzyFit Toolbox

% History:
% 2006/02/09: v1.00, first version.


% if no input argument, use the last fit, which is stored
% in the variable lastfit in the 'base' workspace:
if nargin==0,
    if evalin('base','exist(''lastfit'',''var'')')
        f=evalin('base','lastfit');
    else
        errordlg('No existing fit coefficients. First fit a curve.',...
            'Edit Fit Coefficients','on');
        return;
    end;
end;


if exist('fitparam.m','file'),
    fp=fitparam;
else
    error('No fitparam file found.');
end;


% c is a 2xN cell array containing the parameter
% names and the parameter values
for i=1:length(f.param);
    c{i,1}=f.param{i};
    c{i,2}=f.m(i);       
end;

switch fp.corrcoefmode
    case 'r'
        c{length(f.param)+1,1}='R';
        c{length(f.param)+1,2}=f.r;
    case 'r2',
        c{length(f.param)+1,1}='R^2';
        c{length(f.param)+1,2}=f.r^2;
end;

% displays the coeff in row or in column:
if ~strcmpi(fp.coeffarray,'row'),
    c=c';
end;

% creates a variable called 'FitCoefficients' in the matlab 'base' workspace
% (this is required for using openvar)
assignin('base','FitCoefficients',c)

% calls the Array Editor with this variable "FitCoefficients".
openvar('FitCoefficients');

if nargout==0
    clear c
end;
