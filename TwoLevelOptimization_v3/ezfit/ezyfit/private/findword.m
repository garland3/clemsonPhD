function [posword,word] = findword(inputstr,delim)
%FINDWORD   List of words
%   [POS,WORD] = FINDWORD(INPUTSTR,DELIM) returns all the 'words'
%   (= group of characters separated by delimiters) and their position in
%   the string INPUTSTR. POS is an array, and WORD a cell array of strings.
%
%   F. Moisy
%   Revision: 1.02,  Date: 2012/07/08
%
%   See also SHOWEQBOX, EQ2ML.


% History:
% 2006/02/08: v1.00, first version (was a subfunction of showeqbox)
% 2006/02/17: v1.01, independent m-file.
% 2012/07/08: v1.02, add ',' in the list of delimiters (thanks Patrick)



if nargin==1
    delim='+-*/^()=,'; 
end
rem=inputstr;
nw=0;
while ~isempty(strtok(rem,delim))
    nw=nw+1;
    [word{nw}, rem] = strtok(rem,delim);
    posword(nw) = length(inputstr)-length(rem)-length(word{nw})+1;
end
