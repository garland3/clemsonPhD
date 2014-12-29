function [mass, energy, volume, maxStress] = flywheel(hpoints, VolFractpoints)
% mass, kg
% Joules
% m^3
% Pascals

% Input should be in meters. 


doplot = 1; % Set equal 1 to plot functions


global innerR;
global outerR;
innerR =0.02 % m,  inner shoudl be 0.02
outerR = 0.2 % m,  out should be 0.2

% global p;
% p = 2.78e3;%'density'
global w;
w = 630;%'omega'

% ----------------------
% Model the height with the bezier curve
%----------------------------
deltaRSize = 1/20;



% hpoints = [innerR,0.5;
%             0.04,0.5;
%              0.08,0.5;
%               0.12,0.5;
%              outerR,0.5]
t = linspace(0,1,1/deltaRSize)
h =  bezierInter(hpoints,t)



% ----------------------
% Plot the height and control poiints
%----------------------------
if(doplot ==1)
    subplot(2,2,1);
    plot(h(1,:),h(2,:),'--r')
    hold all
    x =hpoints(:,1);
    y = hpoints(:,2);
    scatter(x,y,150,'rd','fill');
end

% ----------------------
% Model the material volume composition
%----------------------------
%deltaRSize = 1/10; % Resuse, the other one

% Vol Fraction models the Sn in the fly wheel
%global VolFractpoints;
% VolFractpoints = [innerR,0.7;
%             0.04,0.7;
%              0.08,0.7;
%               0.12,0.7;
%              outerR,0.7]
%t = linspace(0,1,1/deltaRSize) % Resuse, the other one
volFrac =  bezierInter(VolFractpoints,t)

% ----------------------
% Plot the volume and control poiints
%----------------------------
%subplot(2,2,1);
if(doplot ==1)
    plot(volFrac(1,:),volFrac(2,:),'--g')
    % hold all
    x =VolFractpoints(:,1);
    y = VolFractpoints(:,2);
    scatter(x,y,150,'gd','fill');
    t= strcat('Green is Vol Fraction, Red is height');
    title(t);

    hold off
end

% ----------------------
% calculate, energy, mass, volume
%----------------------------
volume = CalculateVolume(hpoints, innerR, outerR); % m^3
mass = CalculateMass(hpoints,VolFractpoints, innerR, outerR); % kg
energy = CalculateEnergy(hpoints,VolFractpoints,w, innerR, outerR); % Joule, Newton*meters

% ----------------------
% Calcuale the stresses
%----------------------------

r = h(1,:) % use th interopolated R from the bezier curve.
solution = runEquation(hpoints,VolFractpoints, innerR,outerR);
psi = deval(solution,r);

sigma_r = CalculateRadialStress(psi, h, r);
sigma_theta = CalculateTangentalStress(psi, h,volFrac,r,innerR,outerR);
vonMises = CalculateVonMises(sigma_r,sigma_theta);



maxStress = max(vonMises);

disp(maxStress)

% ----------------------
% Plot Stresses
%----------------------------
if(doplot ==1)
    subplot(2,2,2)
    hold on
    plot(r,sigma_r,'--r')
    plot(r,sigma_theta,'--g')
    plot(r,vonMises,'--b')
    t= strcat('Green=Stress/Sigma Theta: Red=stress radial: Blue is vonMises');
    title(t);

    hold off
end

 
% plot(b(1,:),b(2,:),'--r')
% plot(b(1,:),deriv(2,:),'--g')
% x =hpoints(:,1);
% y = hpoints(:,2);
% scatter(x,y,150,'bd','fill');
% hold off

% psi has 2 rows. h has 2 rows, the second is h, the first is r
function sigma_r = CalculateRadialStress(psi, h, r)
sigma_r = psi(1,:)./(h(2,:).*r);

function sigma_theta = CalculateTangentalStress(psi, h, volFrac, r,~,~)
%global p;
global w;


% t = convertToT(r,innerR,outerR);
% vol = bezierInter(VolFractpoints,t);
p = caculateDensity(volFrac);
p = p(2,:);

sigma_theta = psi(2,:)./h(2,:)   +   p*w^2.*(r.^2);
 

        

function dydx = psi_ODE(r,psi,points,VolFractpoints,innerR,outerR)
v = 0.33; % set poission's ratio to constant for now. 
% global p;
global w;

% Use bezier curve to interpolate the H, and dH needed for this particular
% r value. 
t = convertToT(r,innerR,outerR);
h = bezierInter(points,t);
h = h(2);
dh = DerivBezierInter(points,t);
dh = dh(2);

% Calculate the density, 
% global VolFractpoints;
t = convertToT(r,innerR,outerR);
vol = bezierInter(VolFractpoints,t);
p = caculateDensity(vol);
p = p(2);

dydx = [psi(2);
    (-r*psi(2)+psi(1)-(3+v)*p*w^2*h*r^3+r/h*dh*(r*psi(2)-v*psi(1)))/r^2];

function res = psi_BVs(psi_inner,psi_outer)
res = [psi_inner(1);psi_outer(1)];

function [sol] = runEquation(points,VolFractpoints,innerR,outerR)    
% http://stackoverflow.com/questions/2256229/matlab-how-do-i-pass-a-parameter-to-a-function
solinit = bvpinit(linspace(innerR,outerR,1000),[0 0]);

options = bvpset('RelTol',1e-3);
sol = bvp4c(@(psi,r) psi_ODE(psi,r,points,VolFractpoints,innerR,outerR),@psi_BVs,solinit,options);




%function p = caculateDensity(volFracSn,densitySN,densityAL)


%function p = caculateDensity(volFracSn,densitySN,densityAL)
function v = caculatePoissonsRatio(volFracSn)
% table 3.2 in Huang's disertation, page 46
vSN = 0.3;
vAL =0.3;
v = volFracSn*vSN+(1-volFracSn)*vAL;

function E = caculateYoungsMod(volFracSn)
% table 3.2 in Huang's disertation, page 46
E_SN = 41.4; % GPa
E_AL = 73; % GPa
E = volFracSn*vSN+(1-volFracSn)*vAL;

function vonMises = CalculateVonMises(sigma_r,sigma_theta)

vonMises = (1/2*(sigma_r -sigma_theta).^2+sigma_r.^2+sigma_theta.^2).^(1/2)



