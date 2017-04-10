function remove_efmenu_fig(filenamein,filenameout)
%REMOVE_EFMENU_FIG   Remove the Ezyfit menu from figure files (.FIG)
%   REMOVE_EFMENU_FIG(FILENAMEIN,FILENAMEOUT) removes the Ezyfit menu from
%   the figure file FILENAMEIN, and saves the result under FILENAMEOUT.
%
%   REMOVE_EFMENU_FIG(FILENAMEIN) automatically generates the output
%   filename, by adding the string '_new' to the input filename.
%
%   This function fixes the issue of figure files including the Ezyfit menu
%   opened in a Matlab system running without the Ezyfit toolbox. See
%   EFMENU for details.
%
%   Example:
%     plotsample;   % generate a random plot
%     saveas(gcf,'myfig.fig');
%     remove_efmenu_fig('myfig.fig');
%     % This generates a new file, named 'myfig_new.fig', without the menu
%
%   Acknowledgments to Francis Burton and Nicholas Sinclair, who fixed
%   the issue in the Matlab Central Newsgroup.
%
%   See also EFMENU.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.00,  Date: 2010/07/05
%   This function is part of the EzyFit Toolbox

% History:
% 2010/07/05: v1.00, first version, from nicholas.sinclair@gmail.com

if nargin==1
    [pathstr, name, ext] = fileparts(filenamein);
    filenameout = [name '_new.fig'];
end

v=load(filenamein,'-mat');
hgS_070000=v.hgS_070000;
c=hgS_070000.children;
count = 0;
for i = 1:length(c)
    if isfield(c(i).properties,'Label')
        if(strfind(c(i).properties.Label,'EzyFit'))
            count = count+1;
            ix(count) = i;
        end
    end
end
c(ix)=[];
hgS_070000.children=c;
hgS_070000.properties.CreateFcn = '';
save(filenameout,'hgS_070000');


