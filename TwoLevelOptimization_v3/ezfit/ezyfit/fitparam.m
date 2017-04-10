function fp=fitparam
%FITPARAM  Default settings for the EzyFit toolbox
%   This M-File contains the default settings for the EzyFit Toolbox,
%   such as the extrapolation mode, fit colors, equation box location etc.
%   See the page 'Settings' in the help browser for the list of available
%   settings.
%
%   To change the default settings, edit this file (type: edit fitparam.m)
%   and follow the instructions. You may also choose the item
%   'Default Settings' in the EzyFit menu.
%
%   FP = FITPARAM returns the structure FP containing the settings.
%
%   See also SHOWFIT, EZFIT, SHOWEQBOX, DISPEQFIT, PICKDATA.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.21,  Date: 2007/09/18
%   This function is part of the EzyFit Toolbox




% Fit color mode:
%  - a [R G B] vector specifies a fixed color (eg, [0 0 0] is black)
%  - a string specifies one of eight predefined colors (full name: 'red',
%    'blue' etc, or single letters, 'r', 'b' etc). Type 'doc colorspec' for
%    details.
%  - a numeric value specifies that the color of the fit has the same color
%    as the data, but multiplied by the factor fitcolor (<1 for darker,
%    >1 for lighter)
fp.fitcolor = [1 0.5 0];


% Width and style of the fit (type 'doc linespec' to see the available
% Matlab line specifications):
fp.fitlinewidth = 2;
fp.fitlinestyle = '-';


% Extrapolation mode:
%  'fig'  = extrapolates the fit to the figure limits
%  'data' = extrapolates the fit to the whole data limits, even if only a
%           selection of the data is fitted.
%  'none' = no extrapolation
fp.extrapol = 'none';


% Number of X points used to compute and display the fitted curve:
% (typically 20 to 500, default=200)
fp.npt = 200;


% Display the equation in the command windows: 'on' or 'off'
fp.dispeqmode = 'on';


% Display the equation box in the figure: 'on' or 'off'
fp.dispeqboxmode = 'on';


% Display a legend in the figure: 'on' or 'off'
fp.dispfitlegend = 'off';


% Equation replace mode:
%  'on'    replace each parameter by its numerical value in the equation
%  'off'   keep the parameter names in the equation 
fp.eqreplacemode = 'off';


% Correlation coefficient mode:
%  'r'      display R
%  'r2'     display R^2
%  'none'   display nothing
fp.corrcoefmode = 'r';


% Lin/log mode display: 'on' or 'off'
% (this tells whether Y or LOG(Y) is fitted)
fp.linlogdisp = 'on';


% Open the Array Editor with the fit coefficients after showfit:
% 'on' or 'off'. See Editcoeff for details.
fp.editcoeffmode = 'off';


% Array for the Fit Coefficients: 'row' or 'line':
% (this option is useful for copy-paste the coefficients in Excel)
fp.coeffarray = 'row';


% Automatic call of MakeVarFit after each fit: 'on' or 'off'
% Setting to 'on' creates on the Matlab workspace the variables
% associated to the fit parameters (the variables will be overwritten
% if they already exist!). See Makevarfit for details
fp.automakevarfit = 'off';


% Maximum length of the equation string to be displayed in
% the equation box (longer strings are truncated).
% Set 'maxlength = inf' for no truncation:
fp.maxlengtheq = 35;
%fp.maxlengtheq = inf;


% Location and size of the equation box (in normalized units):
fp.boxlocation = [0.15 0.81 0.3 0.1]; % this is the top-left corner


% Which data to fit when several curves are present in the figure:
% 'first' or 'last'.
fp.whichpickdata = 'first';



% Name of the polynom coefficients. Using a '_' (underscore)
% at the end of the coeff name allows for using the coefficient order
% as a subscript (latex syntax)
% (default: 'a_')
fp.polynom_coeffname = 'a_';


% Number of digits for the coefficient values (default = 5)
fp.numberofdigit = 5;
