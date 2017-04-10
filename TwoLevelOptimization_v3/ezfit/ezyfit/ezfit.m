function f = ezfit(varargin)
%EZFIT   Fit data with arbitrary fitting function
%   EZFIT(FUN) fits the active curve with the function FUN. See below for
%   the syntax of FUN. If FUN is not specified, 'linear' is used.
%
%   By default, the first curve in the active figure is fitted - see
%   FITPARAM to change this default behavior. To fit another curve, select
%   it before calling EZFIT.  If some data are selected by the "Data
%   Brushing" tool (only for Matlab >= 7.6), only those data are fitted.
%
%   EZFIT(X,Y,FUN) or EZFIT(Y,FUN) fits the data (X,Y) (or Y) using the
%   function FUN. X and Y must be vectors of equal length. If X is not
%   specified, X=[1, 2, 3...] is assumed.
%
%   EZFIT(X,Y,FUN), where X is a 1-by-N vector and Y is a 2-by-N matrix,
%   also specifies the weights for Y(2,:) (error bars). By default, when
%   Y is a 1-by-N vector, all the weights are 1.
%
%   Note that EZFIT only computes the coefficients, but does not display the
%   fit. Use SHOWFIT to display the fit.
%
%   The function string FUN can be:
%      - the name of a predefined fitting function (see below).
%      - the name of a user-defined fitting function (see EDITFIT).
%      - an equation, in the form 'y(x)=...', where 'x' represents the
%        X-data, and all the other variables are parameters to be fitted
%        ('a', 'x_0', 'tau', ...). Example: 'y(x)=a*sin(b*x)'. If the
%        left-hand-side 'y(x)' is not specified, 'x' is taken for the
%        X-Data. All the parameter names are accepted, except Matlab
%        reserved strings ('sin', 'pi', ...)
%
%   The predefined fitting functions are:
%      - linear             y = m * x
%      - affine or poly1    y = a*x + b
%      - poly{n}            y = a0 + a1 * x + ... + an * x^n
%      - power              y = c*x^n
%      - sin                y = a * sin (b * x)
%      - cos                y = a * cos (b * x)
%      - exp                y = a * exp (b * x)
%      - log                y = a * log (b * x)
%      - cngauss            y = exp(-x^2/(2*s^2))/(2*pi*s^2)^(1/2)
%      - cfgauss            y = a*exp(-x^2/(2*s^2))
%      - ngauss             y = exp(-(x-x0)^2/(2*s^2))/(2*pi*s^2)^(1/2)
%      - gauss              y = a*exp(-(x-x0)^2/(2*s^2))
%   'ngauss' is a 2-parameters normalized Gaussian, and 'gauss' is a
%   3-parameters non-normalized (free) Gaussian. 'cngauss' and 'cfgauss'
%   are centered normalized and centered free Gaussian, respectively.
%
%   EZFIT is based on Matlab's built-in FMINSEARCH function (Nelder-Mead
%   method), which performs an unconstrained nonlinear minimization of
%   the SSR (sum of squared residuals) with respect to the various
%   parameters. 
%
%   The correlation coefficient R is defined as SSreg / SStot, where
%      SSreg = sum ((y_fit - mean(y)).^2)   % regression sum of squares
%      SStot = sum ((y     - mean(y)).^2)   % total sum of squares
%   (see https://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient)
%   (NB: definition of R changed in Ezyfit v2.44)
%
%   Nonlinear minimization requires starting guesses (or starting estimates)
%   for the fit parameters. By default, all the starting guesses are taken
%   as 1, or, when using predefined fits (e.g., exp, gauss, power...), the
%   starting guesses are determined depending on the range of the data to
%   be fitted. However, in most cases, values closer to the expected
%   result should be specified to "help" the convergence. It is sufficient
%   to choose values that have the correct sign and correct order of
%   magnitude, e.g. 0.01, 1, 100...
%
%   The starting guesses for the parameters may be specified in two ways:
%     - directly in the string FUN, after the fit definition:
%          'c0 + a*sin(pi*x/lambda); c0=1; a=0.1; lambda=100'
%          ('!' or '$' may also be used instead of ';').
%     - by specifying them as an additional input argument for EZFIT:
%          EZFIT(x,y,'c0 + a*sin(pi*x/lambda)',[0.1 1 100]);
%       Note that in this case the parameters must be ordered alphabetically.
%   Note that if both methods are used, only the starting guesses given in
%   the string FUN are considered.
%
%   By default, Y is fitted in linear mode. If you want to fit LOG(Y)
%   instead, you must specify the option 'log' to the string FUN, separated
%   by the symbol ';' or '$' or '!' (eg, FUN='a*x^n;log'). This is
%   specially useful to fit power laws with equally weighted points in a
%   log scale.  If nothing specified, the option 'lin' is used.
%
%   Example:
%        plotsample('power')
%   and compare the results of:
%        ezfit('power;lin')
%        ezfit('power;log')
%
%   F = EZFIT(...) also returns a structure F having the following fields:
%      - name       name of the fit
%      - eq         equation of the fit
%      - param      cell array of strings: names of the parameters
%      - m          values of the coefficients
%      - m0         initial guess for the coefficients
%      - r          correlation coefficient R (Pearson's correlation)
%      - fitmode    'lin' (y is fitted) or 'log' (log(y) is fitted) mode
%
%   This structure F can be further used with SHOWFIT, DISPEQFIT,
%   SHOWEQBOX, MAKEVARFIT and EDITCOEFF.
%
%   From F, you can get the values of the fitted parameters. If you want to
%   create in the current Matlab workspace the variables associated to
%   these parameters, use MAKEVARFIT (or set the option 'automakevarfit'
%   to 'on' in FITPARAM).
%
%   Examples:
%     plotsample damposc
%     f = ezfit('u(t) = c + u_a * sin(2*pi*t/T) * exp(-t/tau); T=5; tau=20');
%     showfit(f);
%     
%     plotsample poly2
%     [x,y] = pickdata;
%     f = ezfit(x, y, 'z(v) = poly3');
%     editcoeff(f);
%
%     plotsample poly2
%     f = ezfit('beta(z) = poly2');
%     showfit(f, 'fitcolor', 'red', 'fitlinewidth', 2);
%
%   See also SHOWFIT, PLOTSAMPLE, DISPEQFIT, EDITCOEFF,
%   FMINSEARCH, MAKEVARFIT, FITPARAM.


%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 2.53,  Date: 2016/04/28
%   This function is part of the EzyFit Toolbox


% History:
% 2005/05/12: v1.00, first version.
% 2005/05/20: v1.10, Use 'eval', with generic functions
% 2005/05/21: v1.11, added the 'lin','log' options.
% 2005/05/24: v1.12, option 'log' by default for 'power' and 'exp'.
% 2005/07/27: v1.13, cosmetics.
% 2005/09/03: v1.14, check arg.
% 2005/10/07: v1.15, gaussian fits added (ngauss and fgauss, centered/not)
% 2005/10/20: v1.16, help text changed.
% 2005/10/31: v1.20, also returns R. cste and poly{n} fits added. Initial
%                    guess defined within the fitting function string. The
%                    order of the output parameters is changed.
% 2005/11/05: v1.21, evaluate strings for initial guess in FUN
% 2005/12/06: v1.22, opens a dialog box if the polynomial order is not
%                    specified.
% 2006/01/13: v1.24, check for negative data in log mode
% 2006/01/19: v1.25, bug fixed from 1.24
% 2006/01/25: v1.26, check the matlab version
% 2006/02/08: v2.00, new syntax. The output argument is now a fit
%                    structure, and the fitting equation string accepts
%                    arbitrary parameter names.
% 2006/02/14: v2.10, lhs 'y(x)=...' now accepted. Now case sensitive
% 2006/03/07: v2.11, bug fixed from v1.25
% 2006/03/09: v2.20, weighted chi-square criterion (ie: error bars accepted
%                    for y) - undocumented
% 2006/09/04: v2.21, '!' and '$' may be used instead of ';' in the FUN
%                    string: this allows to pass the argument like:
%                    fit power!log
% 2006/09/28: v2.22, 'x^1' replaced by 'x' for 1st order polynomial
% 2006/10/18: v2.30, input raw and column vectors accepted. accepts
%                    additional input parameters to change the default settings.
% 2007/04/17: v2.31, help text improved; standard error messages
% 2007/05/16: v2.40, weigthed fits (v2.20) now documented, with bug fixed.
% 2007/08/18: v2.41, help text improved.
% 2007/09/17: v2.50, guess the initial guess m0 for predefined fits
% 2011/11/04: v2.51, bug fixed for fits in log coordinates
% 2012/06/28: v2.52, check version number removed (unstable)
% 2016/04/28: v2.53, definitions of SSE and SSR exchanged (thanks Yoel!!)




% new v1.26, changed v2.11, removed v2.52
% if str2double(version('-release'))<14
%     error('Ezyfit:ezfit:compatibility','EzyFit requires Matlab 7 or higher.');
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit parameters:  (new v2.30)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% loads the default fit parameters:
try
    fp=fitparam;
catch
    error('Ezyfit:ezfit:fitparamNotFound','''fitparam.m'' file not found.');
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

if nargin==0, % if no input argument: do a linear fit of the current curve
    [x,y,h] = pickdata(fp);
    inputstr='linear';
elseif nargin==1
    if ischar(varargin{1}),    %    EZFIT('a*x+b')
        inputstr=varargin{1};
        [x,y,h] = pickdata(fp);
    else                       %    EZFIT(y)
        y=varargin{1};
        x=1:length(y);
        inputstr='linear';
    end
else   % 2 or more input arguments
    if ~isnumeric(varargin{1})           % EZFIT('fun',...)
        inputstr=varargin{1};
        [x,y,h] = pickdata(fp);
        if isnumeric(varargin{2})    % EZFIT('fun',m0,...)
            m0=varargin{2};
        end
    elseif (isnumeric(varargin{1}) && ~isnumeric(varargin{2}))       % EZFIT(y,'a*x+b',...)
        [y,inputstr]=deal(varargin{1:2});
        x=1:length(y);
        if nargin>2
            if isnumeric(varargin{3})   % EZFIT(y,'a*x+b',m0,...)    
                m0=varargin{3};       
            end
        end
    elseif (isnumeric(varargin{1}) && isnumeric(varargin{2}))   % EZFIT(x,y,...)
        [x,y]=deal(varargin{1:2});      
        if nargin>2
            if ischar(varargin{3})
                inputstr=varargin{3};      % EZFIT(x,y,'fun',...)
                if nargin>3
                    if isnumeric(varargin{4})      % EZFIT(x,y,'fun',m0,...)
                        m0=varargin{4};
                    end
                end
            else
                error('Ezyfit:ezfit:syntaxError','Syntax error. 3rd paramater of EZFIT should be a string.');
            end
        end
    end
end               



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some checks about x and y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% turn all input vectors into row (1xN) vectors (new v2.30):
if size(x,1)>size(x,2),
    x=x';
end
if size(y,1)>size(y,2),
    y=y';
end


% check for error bars for y (new v2.20, fixed 2.40)
if size(y,1)>1
    dy = y(2,:);  % second line = error bars (1/weight)
    y = y(1,:);   % first line = data
else
    dy = ones(1,length(y));    % default error bars: 1
end

if length(x)~=length(y),
    error('Ezyfit:ezfit:dimagree','X and Y dimensions must agree.');
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% processing of the string FUN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cleans the input string:
inputstr = strrep(inputstr,' ','');
inputstr = strrep(inputstr,'!',';');
inputstr = strrep(inputstr,'$',';');
inputstr = strrep([inputstr ';'],';;',';'); % ensures that fun is terminated by a ';'.

% the name of the fit is by default the first part of the input string
p=findstr(inputstr,';'); p=p(1);
f.name = inputstr(1:(p-1));

% separates the first part (fitting function itself) and the remainder:
p = strfind(inputstr,';'); p=p(1);
fun = inputstr(1:(p-1));
remfun = inputstr((p+1):end);


% search for predefined fit or user-defined fit
usepredefinedfit = '';
[defaultfit, userfit] = loadfit;
for i=1:length(defaultfit),
    if strcmp(fun, defaultfit(i).name);
        fun = defaultfit(i).eq;
        usepredefinedfit = defaultfit(i).name;   % new v2.50
    end
end
for i=1:length(userfit),
    if strcmp(fun, userfit(i).name),
        fun = userfit(i).eq;
    end
end


% separates again the first part (fitting function itself) and the remainder:
% (because the predefined/user-defined part may itself contain ';')
fun = [strrep(fun,' ','') ';' remfun];
p=strfind(fun,';'); p=p(1);
remfun = fun((p+1):end);
fun = fun(1:(p-1));


% recognize if a lhs is present
peq = strfind(fun,'=');
if ~isempty(peq),
    lhs = fun(1:(peq-1));    % left-hand side
    rhs = fun((peq+1):end);  % right-hand side
else
    lhs = '';
    rhs = fun;
end


% process the lhs (if present)
if ~isempty(lhs),
    pob = strfind(lhs,'('); % position of opening bracket
    pcb = strfind(lhs,')'); % position of closing bracket
    if ~isempty(pob)
        if pob==1,
            f.yvar = 'y';
        else
            f.yvar = lhs(1:(pob-1));
        end
        f.xvar = lhs((pob+1):(pcb-1));
    else
        f.yvar = lhs;
        f.xvar = 'x';
    end
else % if no lhs present:
    f.xvar='x';
    f.yvar='y';
end


% process the 'poly' (polynomial fit) in the rhs
if strfind(rhs,'poly'),   % polynomial fit:
    order=str2double(rhs(5:end));
    if isempty(order), % new v1.22
        str_ord=inputdlg('Order of the polynomial fit','Polynomial order',1,{'2'});
        if ~isempty(str_ord),
            order=str2double(str_ord{1});
            f.name=['poly' str_ord{1}];
        else
            clear f;
            return;
        end
    end
    if order>20,
        error('Ezyfit:ezfit:invalidPolynomDegree','Invalid polynom degree.');
    end
    rhs = [fp.polynom_coeffname '0'];
    for i=1:order,
        if i>1
            rhs = [rhs '+' fp.polynom_coeffname num2str(i) '*' f.xvar '^' num2str(i)];
        else
            rhs = [rhs '+' fp.polynom_coeffname num2str(i) '*' f.xvar];   % new v2.22
        end
    end
end


% search for option 'lin' or 'log':
% (if several are present, use the last one)
% (if none is present, check the Y-scale of the current figure)
% (if no figure present, take 'lin').
fun = strrep([rhs ';' remfun], ';;', ';');
f.fitmode='';
plin=strfind(fun,';lin');
if ~isempty(plin), plin=plin(end); else plin=1; end
plog=strfind(fun,';log');
if ~isempty(plog), plog=plog(end); else plog=1; end
plast=max(plin,plog);
if plast==1, % if nothing specified
    f.fitmode='lin'; % use 'lin'
else
    f.fitmode=fun((plast+1):(plast+3));
end
fun=strrep(fun,';lin','');
fun=strrep(fun,';log','');


% check for negative data in log mode (new v1.23, fixed 1.24):
if (strcmp(f.fitmode,'log') && sum(y<=0)>0)
    disp('Warning: Zero or negative data ignored');
    nonzero=find(y>0);
    x=x(nonzero);
    y=y(nonzero);
end


% separates again the first part (fitting function itself) and the
% remainder:
p = strfind(fun,';'); p=p(1);
f.eq = fun(1:(p-1));
remfun = fun((p+1):end);


 % convert the equation in the matlab syntax (parameters named m(1),
 % m(2)...)
[eqml,param] = eq2ml(f.eq, f.xvar);
maxm = length(param); % number of parameters
if maxm==0  % new v2.20
    error('Ezyfit:ezfit:noParameter','No parameter to be fitted.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% processing the initial guess
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial guess for the m(i) by default (all m(i)=1):
if ~exist('m0','var')
    m0=ones(1,maxm);
    switch usepredefinedfit   % new v2.50 : "magic" initial guesses
        case 'linear'   % m*x
            m0(1) = (y(end)-y(1))/(x(end)-x(1));
        case 'affine'   % a*x+b
            m0(1) = (y(end)-y(1))/(x(end)-x(1));  % a
            m0(2) = mean(y-m0(1)*x);              % b
        case 'affineshift'   % a*(x-b)
            m0(1) = (y(end)-y(1))/(x(end)-x(1));  % a
            m0(2) = mean(x-y/m0(1));              % b
        case 'power'    % a*x^n;
            % if both x and y are of constant sign:
            if sum(diff(sign(x)))==0 || sum(diff(sign(y)))==0
                m0(2) = log(y(end)/y(1)) / log(x(end)/x(1));  % n
                m0(1) = mean(y./(x.^m0(2)));                  % a
            end
        case 'powerc'    % a*x^n+c;
            % if both x and y are of constant sign:
            if sum(diff(sign(x)))==0 || sum(diff(sign(y)))==0
                m0(3) = log(y(end)/y(1)) / log(x(end)/x(1));  % n
                m0(1) = mean(y./(x.^m0(3)));                  % a
            end
            m0(2) = mean(y);  % c
        case 'powershift'    % a*(x+b)^n;
            % if both x and y are of constant sign:
            m0(2) = mean(x);                              % b
            m0(3) = log(y(end)/y(1)) / log((x(end)+m0(2))/(x(1)+m0(2)));  % n
            m0(1) = mean(y./((x+m0(2)).^m0(3)));                  % a           
        case 'exp'    % a*exp(b*x)
            m0(2) = (log(y(end)/y(1))) / (x(end)-x(1));  % b
            m0(1) = mean(y./exp(m0(2)*x));               % a
        case 'expdiv'    % a*exp(x/b)
            m0(2) = (x(end)-x(1)) / (log(y(end)/y(1)));  % b
            m0(1) = mean(y./exp(x/m0(2)));               % a
        case 'explim'    % a*(1-exp(-x/b))
            m0(1) = max(y);                              % a
            m0(2) = mean(-x./(log(1-y/m0(1))));          % b
        case 'expc'    % a*exp(b*x)+c
            m0(3) = mean(y);                             % c
            m0(2) = (log((y(end)-m0(3))/(y(1)-m0(3)))) / (x(end)-x(1));  % b
            m0(1) = mean((y-m0(3))./exp(m0(2)*x));               % a
        case 'log'    % a*log(b*x)
            m0(1) = (y(end)-y(1))/(log(x(end)/x(1)));    % a
            m0(2) = mean(exp(y/m0(1))./x);               % b            
        case 'logc'    % a*log(x)+b
            m0(1) = (y(end)-y(1))/(log(x(end)/x(1)));    % a
            m0(2) = mean(y - m0(1)*log(x));              % b     
        case {'sin','cos'}    % a*sin(b*x), a*cos(b*x)
            m0(1) = std(y,1)*sqrt(2);  % a
            m0(2) = 50/(x(end)-x(1));  % b
        case {'sinc','cosc'}    % a*sin(b*x)+c, a*cos(b*x)+c
            m0(1) = std(y,1)*sqrt(2);  % a
            m0(2) = 50/(x(end)-x(1));  % b
            m0(3) = mean(y);           % c
        case 'sinphi'    % a*sin(b*x+phi)
            m0(1) = std(y,1)*sqrt(2);  % a
            m0(2) = 50/(x(end)-x(1));  % b
            m0(3) = 1;                 % phi
        case 'sinphic'    % a*sin(b*x+phi)+c
            m0(1) = std(y,1)*sqrt(2);  % a
            m0(2) = 50/(x(end)-x(1));  % b
            m0(3) = mean(y);           % c
            m0(4) = 1;                 % phi
        case 'cngauss'
            m0(1) = (mean(x.^2.*y)/mean(y))^(1/2);   % sigma
        case 'cfgauss'   % a*exp(-(x^2)/(2*sigma^2))
            m0(1) = max(y);      % a
            m0(2) = (mean(x.^2.*y)/mean(y))^(1/2);   % sigma
        case 'ngauss'    % exp(-((x-x_0)^2)/(2*sigma^2))/(2*pi*sigma^2)^(1/2)           
            m0(2) = mean(x.*y)/mean(y);             % x_0
            m0(1) = (mean((x-m0(2)).^2.*y)/mean(y))^(1/2);  % sigma
        case {'fgauss','gauss'}    % a*exp(-((x-x_0)^2)/(2*sigma^2))
            m0(1) = max(y);     % a
            m0(3) = mean(x.*y)/mean(y);             % x_0                   
            m0(2) = (mean((x-m0(3)).^2.*y)/mean(y))^(1/2);  % sigma
    end
    % if an initial guess is 0 or infinity or has imaginary part, set it to 1
    for np=1:length(m0)
        if m0(np)==0 || isinf(m0(np)) || imag(m0(np))~=0
            m0(np) = 1;
        end
    end
end
    

% search for initial guess defined into remfun
remfun = strrep([remfun ';'], ';;', ';'); % adds a final ';'
while strfind(remfun,'='),
    peq=strfind(remfun,'='); peq=peq(1);
    pc=strfind(remfun,';'); pc=pc(1);
    if pc>peq,   % if ';' after '='
        par=remfun(1:(peq-1));
        for i=1:maxm,
            if strcmp(param{i},par),
                m0(i)=eval(remfun((peq+1):(pc-1)));
            end
        end
    else
        error('Ezyfit:ezfit:syntaxError','Invalid syntax');
    end
    remfun=remfun((pc+1):end); % removes the processed i.g. and loops back
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fitting itself
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch f.fitmode,
    case 'lin',
        x__ref = x;
        try
            m=fminsearch(@fitlin, m0);
        catch
            error('Ezyfit:ezfit:fminsearchError','Fit: error during the fminsearch procedure');
        end        
        y_fit = eval(eqml);
        % new definition 2016:
        ssreg = sum(abs(y_fit-mean(y)).^2);
        sstot = sum(abs(y-mean(y)).^2);
        f.r = ssreg/sstot;
    case 'log',
        x__ref = x;
        try
          m=fminsearch(@fitlog, m0);   % bug fixed here! (v2.51)
        catch
            error('Ezyfit:ezfit:fminsearchError','Fit: error during the fminsearch procedure');
        end
        y_fit = eval(eqml);
        % new definition 2016:
        ssreg = sum(abs(log(y_fit)-mean(log(y))).^2);
        sstot = sum(abs(log(y)-mean(log(y))).^2);
        f.r = ssreg/sstot;
    otherwise
        error('Ezyfit:ezfit:unknownMode''Unknown fit mode');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% fills the output structure:
f.param = param;
f.m = m;
f.m0 = m0;
f.x = x;
f.y = y;
if sum(dy-1) % store the error bars only if defined
    f.dy = dy;
end
if exist('h','var'), f.hdata=h; end  % handle to the data

% stores the fit in the variable 'lastfit' in the 'base' workspace:
assignin('base','lastfit',f);

if strcmp(fp.automakevarfit,'on')
    makevarfit;
end

% ending displays (if no output argument):
if ~nargout
    if strcmp(fp.dispeqmode,'on') % new v2.30
        dispeqfit(f,fp);
    end
    clear f;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of the main function EZFIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Nested functions that evaluate the fit for prescribed parameters m(i),
% and return the chi2 (sum of the squared difference between the input
% curve and the fit), in lin or log mode:

    function chi2 = fitlin(m)
        y_fit = eval(eqml);
        chi2 = sum(((y_fit - y).^2)./(dy.^2));
    end

    function chi2 = fitlog(m)
        y_fit = eval(eqml);
        chi2 = sum((log(y_fit)-log(y)).^2);
    end
end
