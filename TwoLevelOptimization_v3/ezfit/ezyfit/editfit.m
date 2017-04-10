function editfit(n,name,eq)
%EDITFIT   Edit a user defined fit
%   EDITFIT edit a new user defined fit.
%
%   EDITFIT(N) edits the fit #N.
%
%   EDITFIT(N,NAME,EQ) sets the name NAME and the equation EQ of the user
%   defined fit #N.
%
%   EDITFIT('reset') deletes the user defined fits.
%   EDITFIT('list')  lists the user defined fits
%
%   Example:  editfit(3,'myexp','a*exp(-x/tau); a=0.1; tau=100;');
%
%   See also LOADFIT, EFMENU, EZFIT.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.20,  Date: 2006/02/02
%   This function is part of the EzyFit Toolbox

% History:
% 2006/01/12: v1.00, first version.
% 2006/01/19: v1.10, syntax changed. Now does not ask for the fit number.
% 2006/02/02: v1.20, use a new mat-file for user defined fits


% directory where the ezyfit toolbox is installed:
efroot=fileparts(mfilename('fullpath'));
userfitfile=[efroot filesep 'userfit.mat'];

% option 'reset' or 'list':
if (nargin==1) && (~isnumeric(n))
    if strcmp(n,'reset'),
        clear userfit;
        if exist(userfitfile,'file'),
            delete(userfitfile);
        end;
        loadfit('user');
        return;
    elseif strcmp(n,'list'),
        userfit=loadfit('user');
        for i=1:length(userfit),
            disp(['Fit #' num2str(i) ': ' userfit(i).name]);
            disp(['        ' userfit(i).eq]);
        end;
        return;
    end;
end;

userfit=loadfit('user');

% 3 arguments: sets the fit #n to its new definition.
if (nargin==3),
    if n<=(length(userfit)+1),
        userfit(n).name=name;
        userfit(n).eq=eq;
        save(userfitfile,'userfit');
        efmenu;  % refresh the ezyfit menu
        return;
    else
        error('Wrong fit number.');
    end;
end;


% no argument: creates a new fit.
if (nargin==0)
    n=length(userfit)+1;
    userfit(n).name = ['fit' num2str(n)];
    userfit(n).eq = 'a*x+b*x^2 ; a=1 ; b=0.1';
end;

if n>length(userfit),
    error('Wrong fit number.');
end;
    
answer = inputdlg({'Name (enter an empty string to delete this fit):',...
    'Equation y(x) = a*x + b*x^2...; a=1; ...'},...
    ['Edit User Fit #' num2str(n)], 1,...
    {userfit(n).name,userfit(n).eq});

if ~isempty(answer)
    if isempty(answer{1}), % deletes this fit
        userfit = userfit([1:(n-1) (n+1):end]);
    else
        userfit(n).name=answer{1};
        userfit(n).eq=answer{2};
    end;
    save(userfitfile,'userfit');
    efmenu;  % refresh the ezyfit menu
else
    return;
end;

