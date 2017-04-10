function [fit1, fit2] = loadfit(opt)
%LOADFIT   Load the predefined and the user-defined fitting functions.
%   [DEFAULTFIT, USERFIT] = LOADFIT loads the predefined and the user-
%   defined fitting functions.
%
%   DEFAULTFIT = LOADFIT('default') loads only the predefined fits.
%   USERFIT = LOADFIT('user') loads only the user-defined fits.
%
%   DEFAULTFIT and USERFIT are structure arrays containing two fields,
%   'name' and 'eq'.
%
%   If the file(s) for the predefined and/or the user-defined fits do(es)
%   not exist, LOADFIT creates it (them).
%
%   See also EDITFIT, EFMENU.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.10,  Date: 2007/09/18
%   This function is part of the EzyFit Toolbox

% History:
% 2006/02/02: v1.00, first version.
% 2007/09/18: v1.10, new default fits (expc, sinc etc)

if nargin==0, opt='defaultuser'; end

user=0;
default=0;

% directory where the ezyfit toolbox is installed:
efroot=fileparts(mfilename('fullpath'));

if strfind(lower(opt),'default'),
    default=1;
    % load the predefined fits:
    defaultfitfile=[efroot filesep 'defaultfit.mat'];
    if exist(defaultfitfile,'file'),
        load(defaultfitfile);
    else
        defaultfit(1).name = 'linear';          defaultfit(1).eq = 'a*x';
        defaultfit(2).name = 'affine';          defaultfit(2).eq = 'a*x+b';
        defaultfit(3).name = 'affineshift';     defaultfit(3).eq = 'a*(x-b)';
        defaultfit(4).name = 'power';           defaultfit(4).eq = 'a*x^n';
        defaultfit(5).name = 'powerc';          defaultfit(5).eq = 'a*x^n+c';
        defaultfit(6).name = 'exp';             defaultfit(6).eq = 'a*exp(b*x)';
        defaultfit(7).name = 'expc';            defaultfit(7).eq = 'a*exp(b*x)+c';
        defaultfit(8).name = 'log';             defaultfit(8).eq = 'a*log(b*x)';
        defaultfit(9).name = 'logc';            defaultfit(9).eq = 'a*log(x)+b';
        defaultfit(10).name = 'sin';            defaultfit(10).eq = 'a*sin(b*x)';
        defaultfit(11).name = 'sinc';           defaultfit(11).eq = 'a*sin(b*x)+c';
        defaultfit(12).name = 'cos';            defaultfit(12).eq = 'a*cos(b*x)';
        defaultfit(13).name = 'cosc';           defaultfit(13).eq = 'a*cos(b*x)+c';
        defaultfit(14).name = 'sinphi';         defaultfit(14).eq = 'a*sin(b*x+phi)';
        defaultfit(15).name = 'sinphic';        defaultfit(15).eq = 'a*sin(b*x+phi)+c';
        defaultfit(16).name = 'cngauss';        defaultfit(16).eq = 'exp(-(x^2)/(2*sigma^2))/(2*pi*sigma^2)^(1/2)';
        defaultfit(17).name = 'cfgauss';        defaultfit(17).eq = 'a*exp(-(x^2)/(2*sigma^2))';
        defaultfit(18).name = 'ngauss';         defaultfit(18).eq = 'exp(-((x-x_0)^2)/(2*sigma^2))/(2*pi*sigma^2)^(1/2)';
        defaultfit(19).name = 'fgauss';         defaultfit(19).eq = 'a*exp(-((x-x_0)^2)/(2*sigma^2))';
        defaultfit(20).name = 'gauss';          defaultfit(20).eq = 'a*exp(-((x-x_0)^2)/(2*sigma^2))';
        defaultfit(21).name = 'expdiv';         defaultfit(21).eq = 'a*exp(x/b)';
        defaultfit(22).name = 'explim';         defaultfit(22).eq = 'a*(1-exp(-x/b))';
        defaultfit(23).name = 'powershift';     defaultfit(23).eq = 'a*(x+b)^n';
        
        save(defaultfitfile,'defaultfit');
    end
end


if strfind(lower(opt),'user')
    user=1;
    % load the user-defined fits:
    userfitfile=[efroot filesep 'userfit.mat'];
    if exist(userfitfile,'file'),
        load(userfitfile);
    else
        userfit(1).name = 'fit1';       userfit(1).eq = 'a*exp(-x/tau); a=1; tau=10';
        userfit(2).name = 'fit2';       userfit(2).eq = 'c+(x/x_0)^n';
        userfit(3).name = 'fit3';       userfit(3).eq = 'E(k)=C*k^(-n)*exp(-k/k_c); log; C=5; n=2; k_c=100';

        save(userfitfile,'userfit');
    end
end



% output arguments:
if user && ~default,
    fit1=userfit;
    fit2=[];
elseif default && ~user,
    fit1=defaultfit;
    fit2=[];
elseif default && user,
    fit1=defaultfit;
    fit2=userfit;
end
