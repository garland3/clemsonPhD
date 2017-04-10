function str=greekize(inputstr)
%GREEKIZE   Greekize a string
%
%   F. Moisy
%   Revision: 1.01,  Date: 2006/02/17
%
%   See also SHOWEQBOX.


% History:
% 2006/02/08: v1.00, first version (was a subfunction of showeqbox)
% 2006/02/17: v1.01, independent m-file.


greek={'alpha','beta','gamma','delta','epsilon','zeta','eta',...
    'theta','vartheta','iota','kappa','lambda','mu','nu','xi','pi','varpi',...
    'rho','sigma','varsigma','tau','upsilon','phi','chi','psi','omega',...
    'Gamma','Delta','Theta','Lambda','Xi','Pi','Sigma','Upsilon','Phi',...
    'Psi','Omega'};

str=strrep(inputstr,' ','');
[posword,word] = findword(str,'+-*/^()=_'); % ('_' added, to interpret subscripted greek letters)
for n=1:length(greek),
    if ~isempty(strmatch(greek{n},word)),
        % greekize only if the greek letter is a parameter
        % (ie, do NOT greekize 'peta' in 'p\eta')
        str=strrep(str,greek{n},['\' greek{n}]);
    end;
end;

% various latex replacements:
str=strrep(str,'*',' '); % use '\cdot' or '\times' or ' ' here
str=strrep(str,'/',' / ');
str=strrep(str,'+-','-');
str=strrep(str,'+',' + ');
str=strrep(str,'-',' - ');
str=strrep(str,'^','\^');
str=strrep(str,'=',' = ');
str=strrep(str,'  ',' ');
