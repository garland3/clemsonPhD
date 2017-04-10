function about_ef(opt)
%ABOUT_EF display the "About" information of the EzyFit toolbox
%    ABOUT_EF displays the dialog box 'About EzyFit'.
%    ABOUT_EF('command') displays the 'about' info in the command window.
%
%   See also EFMENU, CHECKUPDATE_EF.


%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.02,  Date: 2006/10/31
%   This function is part of the EzyFit Toolbox


% History:
% 2005/12/12: v1.00, first version.
% 2005/04/10: v1.01, Toolbox renamed 'EzyFit'
% 2006/10/31: v1.02, modal window

% gr_dummy not defined for old versions
if verLessThan('matlab','8.4')
    eval('gr_dummy = 0;');
else
    gr_dummy = groot;
end

if nargin==0
    opt='command';
end

v=ver('ezyfit');

switch opt
    case 'command'
        disp('EzyFit');
        disp('A free curve fitting toolbox for Matlab');
        disp(['Version ' v.Version ' (' v.Date ')']);
        disp('Frederic Moisy');
    otherwise
        a=imread('about_ef.jpg');
        ss=get(gr_dummy,'ScreenSize');
        figure('Position',[(ss(3)-size(a,2))/2 (ss(4)-size(a,1))/2 size(a,2) size(a,1)]);
        image(a); set(gca,'Position',[0 0 1 1]);
        axis off;
        set(gcf,'Toolbar','none');
        set(gcf,'Menubar','none');
        set(gcf,'Numbertitle','off');
        set(gcf,'Name','About EzyFit');
        set(gcf,'Resize','off');
        set(gcf,'WindowStyle','modal');
        delete(findobj(gcf,'Label','Ezyfit'));
        
        annotation('textbox',[.175 .33 .4 .2],'Color',.65*[1 1 1],'LineStyle','none',...
            'FontName','Verdana','FontSize',8,'FontWeight','bold',...
            'String',['Version ' v.Version ' (' v.Date ')']);
end
