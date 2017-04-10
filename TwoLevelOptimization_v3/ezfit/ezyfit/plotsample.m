function [x,y]=plotsample(opt,varargin)
%PLOTSAMPLE  Display a sample plot.
%   PLOTSAMPLE displays some noisy sample data. Try to fit the data with
%   SHOWFIT or SELECTFIT, following the instructions given in the command
%   window. The sample data is chosen randomly among the 10 predefined
%   plots described below.
%
%   PLOTSAMPLE is also available from the item 'Plot Sample' of the EzyFit
%   menu (see EFMENU).
%
%   PLOTSAMPLE(OPT)  specifies the sample plot:
%     'power':   noisy power law
%     'linear':  noisy data to be fitted by a linear function
%     'osc':     oscillations
%     'damposc': damped oscillations
%     'cste':    noisy constant
%     'exp':     noisy exponential decay
%     'hist':    histogram of 1000 realizations of a random variable
%     'hist2':   histogram with two gaussian peaks.
%     'powco':   noisy power law with an exponential cut-off.
%     'poly2':   3 curves to be fitted by a 2nd order polynomial fit
%
%   PLOTSAMPLE(N) specifies the number of the sample plot (1 to 10).
%
%   PLOTSAMPLE(...,'fit') also displays the fit.
%   PLOTSAMPLE(...,'nodisp') does not display the help text in the command
%   window.
%
%   [X,Y] = PLOTSAMPLE(...) also returns the sample data (this is
%   equivalent to PLOTSAMPLE; [X,Y] = PICKDATA;).
%
%   Example:
%       plotsample power fit
%       (it shows a power law and fits it).
%
%   See also SHOWFIT, EZFIT, PICKDATA, EFMENU.

%   F. Moisy, moisy_at_fast.u-psud.fr
%   Revision: 1.36,  Date: 2009/09/12
%   This function is part of the EzyFit Toolbox

% History:
% 2005/10/17: v1.00, first version.
% 2005/10/31: v1.10, option added. output argument (x,y).
% 2005/11/02: v1.11, another plot sample.
% 2005/12/06: v1.20, displays the Ezyfit menu.
% 2005/12/13: v1.21, displays the text with a matlab link
% 2006/02/03: v1.30, fits with free parameter names; new plot samples.
% 2006/02/10: v1.31, also clears the annotation objects
% 2006/02/16: v1.32, fits with free lhs 'y(x)=...'. Simpler plot title.
%                    option fitmode added.
% 2006/03/08: v1.33, the figure is not reset any more; this caused bugs
%                    with insert line.
% 2006/09/05: v1.34, bug in the help text fixed.
% 2006/10/29: v1.35, 'hist' now plots a true histogram plot
% 2007/09/12: v1.36, new option 'nodisp' and minor changes

% error(nargchk(0,2,nargin));

rmfit; % remove previous fits

listopt = {'power','linear','osc','damposc','cste','exp','hist','hist2','powco','poly2'};
nsample = length(listopt);

if nargin==0
    opt=floor(nsample*rand+1);
end

if isnumeric(opt)
    opt=listopt{max(1,min(opt,nsample))};
end

switch lower(opt)
    case 'power'
        np = 10+round(50*rand);
        x=logspace(1,3,np)*exp((rand-0.5)*6);
        x=(1 + 0.3*(rand(1,np)-0.5)).*x;
        expo = (0.5+rand*2)*sign(rand-0.5);
        y=exp((rand-0.5)*6)*rand*(1 + .5*(rand(1,np)-0.5)).*x.^expo;
        loglog(x, y, 'bo');
        axisl
        xlabel('x'); ylabel('y');
        f='power';
        txtl=['Try to fit with <a href="matlab:showfit(''' f ''')">showfit(''' f ''')</a>'];
    case 'linear'
        x=linspace(0,12,20);
        y=(rand*1.4+0.9)*x+(rand-0.5)*20;
        x=x+(rand(1,20)-.5)*1.2;
        y=y+(rand(1,20)-.5)*1.6;
        plot(x,y,'ks');
        xlabel('x'); ylabel('y');
        f='affine';
        txtl='Try to fit with <a href="matlab:showfit(''poly1'')">showfit(''poly1'')</a> or <a href="matlab:showfit(''affine'')">showfit(''affine'')</a>';
    case 'osc'
        np = 200 + round(600*rand);
        x=linspace(0,40,np);
        y=2*rand+(2*rand+1)*(1+0.2*rand(1,np)).*sin(rand + x/(1+.1*rand));
        plot(x,y,'b.');
        xlabel('t (seconds)'); ylabel('U(t) [Volts]');
        f='U(t)=U_0+A*sin(omega*t + phi); U_0=1; omega=1; A=2; phi=0.5';
        txtl=['Try to fit with <a href="matlab:showfit(''' f ''')">showfit(''' f ''')</a>'];     
    case 'damposc'
        np = 200 + round(600*rand);
        x=linspace(0,40,np);
        y=2*rand+(2*rand+1)*(1+0.7*rand(1,np)).*sin(x/(1+.3*rand)).*exp(-x/(8+.3*rand));
        plot(x,y,'b.');
        xlabel('t (seconds)'); ylabel('U(t) [Volts]');
        f='U(t)=offset+U_0*sin(2*pi*t/T)*exp(-t/tau_0); offset=1; T=5; U_0=2; tau_0=10';
        txtl=['Try to fit with <a href="matlab:showfit(''' f ''')">showfit(''' f ''')</a>'];
    case 'cste'
        x=linspace(0,20,20);
        y=(rand+.5)*11.1+(rand(1,20)-.5)*(rand+1)*1.3;
        x=x+(rand(1,20)-.5)*.7;
        plot(x,y,'ro');
        axis([0 20 0 1.3*max(y)]);
        xlabel('x'); ylabel('y');
        f='cste';
        txtl='Try to fit with <a href="matlab:showfit(''cste'')">showfit(''cste'')</a> or <a href="matlab:showfit(''affine'')">showfit(''affine'')</a>';
    case 'exp'
        np=15+round(rand*30);
        x=linspace(0,10+round(20*rand),np);
        y=(rand+0.5)*3*exp(-x/(3.8*rand+2));
        x=abs(x+(rand(1,np)-.5)*.7);
        y=abs(y+(rand(1,np)-.5)*0.05);
        plot(x,y,'k*');
        xlabel('t (s)'); ylabel('N(t)');
        f='N(t)=N_0*exp(-t/tau)';
        txtl=['Try to fit with <a href="matlab:showfit(''' f ''')">showfit(''' f ''')</a>'];
    case 'hist'
        np = 300 + round(3000*rand);
        a=ones(1,np);
        for i=1:np
            a(i)=rand+rand+rand+rand+rand+rand+rand+rand-4;
        end;
        a=a+(rand-0.5)*2;  % shift the center
        x=linspace(-5,5,15+round(rand*100));
        scal=exp(8*(rand-0.5));
        hist(a*scal,x*scal);
        xlabel('x'); ylabel('hist');
        f='gauss';
        txtl=['Try to fit with <a href="matlab:showfit(''' f ''')">showfit(''' f ''')</a>'];
    case 'hist2'
        np = 7000 + round(rand*3000);
        a=ones(1,np);
        for i=1:round(np/3),
            a(i)=(rand+rand+rand+rand+rand+rand+rand+rand)/2+5;
        end;
        for i=(round(np/3)+1):np,
            a(i)=(rand+rand+rand+rand+rand+rand+rand+rand)*1.7+8;
        end;
        a((np+1):(np+500))=rand(1,500)*20;
        x=linspace(0,20,400+round(rand*100));
        y=hist(a,x);
        plot(x,y,'r+');
        xlabel('x'); ylabel('hist');
        f='a_1*exp(-(x-m_1)^2/(2*s_1^2))+a_2*exp(-(x-m_2)^2/(2*s_2^2)); a_1=100; a_2=100; m_1=8; m_2=16';
        txtl=['Try to fit with <a href="matlab:showfit(''' f ''')">showfit(''' f ''')</a>'];
    case 'powco'
        x=logspace(0,3,50);
        x=(1 + 0.4*(rand(1,50)-0.5)).*x;
        y=4*rand*(1 + .9*(rand(1,50)-0.5)).*x.^(-(rand+1)).*exp(-x/(100*(rand+.5)));
        loglog(x, y, 'd');
        xlabel('k'); ylabel('E(k)');
        f='E(k)=C*k^(-n)*exp(-k/k_c); log; C=5; n=2; k_c=100';
        txtl=['Try to fit with <a href="matlab:showfit(''' f ''')">showfit(''' f ''')</a>'];
    case 'poly2',
        x=linspace(0,15,20);
        y=rand*2.2*x+(rand*0.6-0.35)*x.^2+(rand(1,20)-.5)*2.5;
        y2=rand*2.2*x+(rand*0.6-0.35)*x.^2+(rand(1,20)-.5)*2.5;
        y3=rand*2.2*x+(rand*0.6-0.35)*x.^2+(rand(1,20)-.5)*2.5;
        plot(x,y,'bo',x,y2,'rs',x,y3,'md');
        xlabel('x'); ylabel('y');
        f='poly2';
        txtl=['Click on a curve, and try to fit with <a href="matlab:showfit(''' f ''')">showfit(''' f ''')</a>'];
    otherwise
        error('Invalid plotsample option.');
end;
title('Sample data');
efmenu;

if any(strcmpi(varargin,'fit'))
    showfit(f);    
else
    if ~any(strncmpi(varargin,'nodisp',2))  % new v1.36
        disp(txtl);
    end
end

if nargout==0
    clear x y
end
