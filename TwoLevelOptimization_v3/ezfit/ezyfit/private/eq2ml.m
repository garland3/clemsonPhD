function [eqml,param]=eq2ml(eq,xvar)
%EQ2ML Convert a free parameter equation to a Matlab equation
%   [EQML, PARAM] = EQ2ML(EQ,XVAR) converts the free-parameter equation
%   string EQ into a Matlab equation string for further evaluation with
%   EVAL or FMINSEARCH, in which the free parameters are replaced by 'm(1)',
%   'm(2)' (in alphabetical order) and the 'x' variable XVAR is replaced by
%   'x__ref'. The names of the parameters are returned in the cell array
%   PARAM. If XVAR is not specified, 'x' is taken by default. In addition,
%   the equation string is 'vectorized' (i.e. '*', '/' and '^' are replaced
%   by '.*', './', '.^').
%
%   Example:  [eqml, param] = eq2ml('c+t_0/t','t')
%      returns 'm(1)+m(2)./x__ref' and {'c','t_0'}
%
%   See also EZFIT, SELECTFIT, SHOWFIT.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.11,  Date: 2006/03/10
%   This function is part of the EzyFit Toolbox

% History:
% 2006/02/06: v1.00, first version.
% 2006/02/13: v1.10, free 'x' variable name.
% 2006/03/10: v1.11, an empty cell array is returned if no param present.


if nargin==1,
    xvar='x';
end

% equation in the matlab syntax (parameters named m(1), m(2)...)
eqml = strrep(eq,' ',''); 

% param0 contains all variables, including the 'x' variable;
% param contains only the parameters of the fit, excluding the 'x'
% variable.

param0=argnames(inline(eqml));
param={};
maxm=0;
for i=1:length(param0),
    if ~(strcmp(param0{i},xvar)),       
        maxm = maxm+1;
        param{maxm} = param0{i};
    end
end


% makes the replacement:
% (before replacing a parameter by 'm(i)', first check
% that the token actually consists in a group of characters
% separated by delimiters. the list of token has to be
% re-computed after each replacement).

[posword,word] = findword(eqml);
for i=1:maxm,
    for j=length(word):-1:1,
        if strcmp(word{j}, param{i}),
            eqml = [eqml(1:(posword(j)-1)) '#(' num2str(i) ')' eqml((posword(j)+length(param{i})):end)];
        elseif strcmp(word{j}, xvar),
            eqml = [eqml(1:(posword(j)-1)) 'x__ref' eqml((posword(j)+length(xvar)):end)];
        end
    end
    [posword,word] = findword(eqml);
end
eqml=strrep(eqml,'#','m'); % replace by '#(xx)', and after by 'm(xx)', so that
                           % the letter 'm' can be used as a parameter
                           % name too.

% inserts a '.' before any '^', '*' or '/' 
eqml=vectorize(eqml);

